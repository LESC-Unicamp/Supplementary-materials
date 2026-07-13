! ************************************************************************************************ !
!                                       PRESSURE CALCULATION                                       !
! ************************************************************************************************ !
!             This subroutine is used to calculate the pressure and the derivatives of             !
!                      the Helmholtz free energy with respect to the volume.                       !
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
! Calculates pressure and the derivatives of the Helmholtz free energy with respect to the volume  !
! for mixtures                                                                                     !
! ************************************************************************************************ !
SUBROUTINE Calculate_Pressure( mFraction, mVolume, Temperature, Pressure, IsothermalCompressibility, ThermalExpansionCoefficient, &
&                              CompressibilityFactor )

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
REAL( Kind= Real64 )                                           :: Pressure                                   ! Pressure
REAL( Kind= Real64 )                                           :: IsothermalCompressibility                  ! Isothermal compressibility
REAL( Kind= Real64 )                                           :: ThermalExpansionCoefficient                ! Thermal expansion coefficient
REAL( Kind= Real64 )                                           :: CompressibilityFactor                      ! Compressibility factor
REAL( Kind= Real64 )                                           :: mVolume                                    ! Molar volume
REAL( Kind= Real64 )                                           :: Temperature                                ! Temperature
REAL( Kind= Real64 )                                           :: AngleAverageMixture                        ! Angle average of the mixture
REAL( Kind= Real64 )                                           :: HSSecondVirialCoefficientMixture           ! Second virial coefficient of a mixture of hard spheres
REAL( Kind= Real64 )                                           :: SecondVirialCoefficientMixtureSingle       ! Second virial coefficient of a mixture of non-spherical rigid bodies
REAL( Kind= Real64 )                                           :: SecondVirialCoefficientMixture             ! Second virial coefficient of a mixture of non-spherical rigid bodies
REAL( Kind= Real64 )                                           :: rNumberDensity                             ! Reduced number density
REAL( Kind= Real64 )                                           :: dNumberDensity_dVolume                     ! First derivative of the number density with respect to the volume
REAL( Kind= Real64 )                                           :: dInverseNumberDensity_dVolume              ! First derivative of the inverse number density with respect to the volume
REAL( Kind= Real64 )                                           :: d2NumberDensity_d2Volume                   ! Second derivative of the number density with respect to the volume
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
REAL( Kind= Real64 )                                           :: dIdealFreeEnergy_dVolume                   ! First derivative of the ideal contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dpFirstOrderCoefficient_dVolume            ! First derivative of the first-order perturbation coefficient with respect to the volume
REAL( Kind= Real64 )                                           :: dpSecondOrderCoefficient_dVolume           ! First derivative of the second-order perturbation coefficient with respect to the volume
REAL( Kind= Real64 )                                           :: dHSBoublikMansoori_dVolume                 ! First derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume
REAL( Kind= Real64 )                                           :: dHCBBoublik_dVolume                        ! First derivative of the Boublik hard convex body term with respect to the volume
REAL( Kind= Real64 )                                           :: dExcludedVolumeFreeEnergy_dVolume          ! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dHSIsothermalCompressibility_dVolume       ! First derivative of the hard-sphere isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                                           :: dHCBIsothermalCompressibility_dVolume      ! First derivative of the hard-convex body isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                                           :: dMeanAttractiveFreeEnergy_dVolume          ! First derivative of the mean-attractive energy contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dMeanFluctuationFreeEnergy_dVolume         ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dMonomerFreeEnergy_dVolume                 ! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dMeanFluctuationFreeEnergy_dTemperature    ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dTotalFreeEnergy_dVolume                   ! First derivative of the total Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2IdealFreeEnergy_d2Volume                 ! Second derivative of the ideal contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2ExcludedVolumeFreeEnergy_d2Volume        ! Second derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2HSBoublikMansoori_d2Volume               ! Second derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume
REAL( Kind= Real64 )                                           :: d2HCBBoublik_d2Volume                      ! Second derivative of the Boublik hard convex body term with respect to the volume
REAL( Kind= Real64 )                                           :: d2MeanAttractiveFreeEnergy_d2Volume        ! Second derivative of the mean-attractive energy contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2pFirstOrderCoefficient_d2Volume          ! Second derivative of the first-order perturbation coefficient with respect to the volume
REAL( Kind= Real64 )                                           :: d2HSIsothermalCompressibility_d2Volume     ! Second derivative of the hard-sphere isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                                           :: d2HCBIsothermalCompressibility_d2Volume    ! Second derivative of the hard-convex body isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                                           :: d2MeanFluctuationFreeEnergy_d2Volume       ! Second derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2pSecondOrderCoefficient_d2Volume         ! Second derivative of the second-order perturbation coefficient with respect to the volume
REAL( Kind= Real64 )                                           :: d2MonomerFreeEnergy_d2Volume               ! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2TotalFreeEnergy_d2Volume                 ! Second derivative of the total Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dxIdealFreeEnergy_dVolume_dTemperature     ! Cross derivative of the ideal contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dxEVFreeEnergy_dVolume_dTemperature        ! Cross derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dxSecondOrderTPT_dVolume_dTemperature      ! Cross derivative of the second-order perturbation coefficient with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dxMonomerFEnergy_dVolume_dTemperature      ! Cross derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dxTotalFreeEnergy_dVolume_dTemperature     ! Cross derivative of the total Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dHigherOrderTerms_dVolume                  ! First derivative of the higher-order terms with respect to the volume
REAL( Kind= Real64 )                                           :: d2HigherOrderTerms_d2Volume                ! Second derivative of the higher-order terms with respect to the volume
REAL( Kind= Real64 )                                           :: dHigherOrderFEnergy_dVolume                ! First derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: d2HigherOrderFEnergy_d2Volume              ! Second derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                                           :: dxHigherOrderTerms_dVolume_dTemperature    ! Cross derivative of the higher-order terms with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: dxHigherOrderFEnergy_dVolume_dTemperature  ! Cross derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                                           :: HigherOrderTerms                           ! Higher-order perturbation terms
REAL( Kind= Real64 )                                           :: Factorial                                  ! Factorial function
REAL( Kind= Real64 )                                           :: mHigherOrderHelmholtzFreeEnergy            ! Higher-order contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: AspectRatioMixture                         ! Aspect ratio of the mixture
REAL( Kind= Real64 )                                           :: sumVolume, sumRadius, sumArea, sumSqRadius ! Morphological descriptors (sum)
REAL( Kind= Real64 )                                           :: sumVolumeSurface                           ! Morphological descriptors (sum)
REAL( Kind= Real64 )                                           :: geoAux, dAux1, dAux2                       ! Auxiliars
REAL( Kind= Real64 )                                           :: HCBAux1, HCBAux2, HCBAux3                  ! Auxiliars
REAL( Kind= Real64 )                                           :: mNonSphericityMixture                      ! Non-sphericity parameter of the mixture
REAL( Kind= Real64 )                                           :: DiameterSphereMixture                      ! Diameter of the representative sphere of the mixture
REAL( Kind= Real64 )                                           :: DiameterSphereMixtureCubic                 ! Diameter of the representative sphere of the mixture (cubic)
REAL( Kind= Real64 )                                           :: ZhangFactor                                ! Zhang's correction factor
REAL( Kind= Real64 )                                           :: gAux1, gAux2                               ! Auxiliars (nonsphericity expressions)
REAL( Kind= Real64 )                                           :: dZhangFactor_dVolume                       ! First derivative of the Zhang's correction factor with respect to the volume
REAL( Kind= Real64 )                                           :: d2ZhangFactor_d2Volume                     ! Second derivative of the Zhang's correction factor with respect to the volume
REAL( Kind= Real64 )                                           :: rDensityMixture                            ! Reduced density of the mixture
REAL( Kind= Real64 ), DIMENSION( 3 )                           :: dAuxiliaryHSFactor_dVolume                 ! Auxiliary factor of a first derivative with respect to the volume (hard-sphere contribution)
REAL( Kind= Real64 ), DIMENSION( 3 )                           :: d2AuxiliaryHSFactor_d2Volume               ! Auxiliary factor of a second derivative with respect to the volume (hard-sphere contribution)
REAL( Kind= Real64 ), DIMENSION( 6 )                           :: d2AuxHSIC_d2Volume                         ! Auxiliary factor of a second derivative with respect to the volume (hard-sphere isothermal compressibility)
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: drDensity_dVolume                          ! First derivative of the reduced densities with respect to the volume
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: d2rDensity_d2Volume                        ! Second derivative of the reduced densities with respect to the volume
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: rDensity                                   ! Reduced densities
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: rDensityConstants                          ! Reduced density constants
REAL( Kind= Real64 ), DIMENSION( 2, 9 )                        :: CSWAlphaCoefficients                       ! Coefficients of the nonsphericity (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3, 9 )                        :: CSW1CoefficientsMatrix                     ! Matrix of coefficients of the effective packing fraction (first-order coefficient)
REAL( Kind= Real64 ), DIMENSION( 3, 9 )                        :: CSW2CoefficientsMatrix                     ! Matrix of coefficients of the effective packing fraction (second-order coefficient)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: MomentInertia                              ! Moment of inertia (non-spherical rigid bodies)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Translational         ! Thermal de Broglie wavelength (translational contribution)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Translational_Cb      ! Thermal de Broglie wavelength (translational contribution, cubic)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Rotational            ! Thermal de Broglie wavelength (rotational contribution)
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
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFraction_dVolume                ! First derivative of the effective packing fraction with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionYK_dVolume              ! First derivative of the effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionSL_dVolume              ! First derivative of the effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionCSW_dVolume             ! First derivative of the effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffPFraction_dVolume_dDensity            ! Cross derivative of the effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffPFractionYK_dVolume_dDensity          ! Cross derivative of the effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffPFractionSL_dVolume_dDensity          ! Cross derivative of the effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffPFractionCSW_dVolume_dDensity         ! Cross derivative of the effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffectiveRDF_dVolume_dDensity            ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffectiveRDFYukawa_dVolume_dDensity      ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffectiveRDFSutherland_dVolume_dDensity  ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxEffectiveRDFCSW_dVolume_dDensity         ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcEffectiveRDF_dVolume                     ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcEffectiveRDFYukawa_dVolume               ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcEffectiveRDFSutherland_dVolume           ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcEffectiveRDFCSW_dVolume                  ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcMeanAttractiveEnergy_dVolume             ! First derivative of the mean-attractive energy with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcMeanAttractiveEnergySutherland_dVolume   ! First derivative of the mean-attractive energy with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcMeanAttractiveEnergyYukawa_dVolume       ! First derivative of the mean-attractive energy with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcMeanAttractiveEnergyCSW_dVolume          ! First derivative of the mean-attractive energy with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcMeanEnergyFluctuations_dVolume           ! First derivative of the mean-attractive energy fluctuations with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dxFirstTPTCoefficient_dVolume_dDensity     ! Cross derivative of the first-order perturbation coefficient with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dcRDFunction_dVolume                       ! First derivative of the contact radial distribution function with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cMeanAttractiveEnergy_d2Volume           ! Second derivative of the mean-attractive energy with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cMeanAttractiveEnergySL_d2Volume         ! Second derivative of the mean-attractive energy with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cMeanAttractiveEnergyYK_d2Volume         ! Second derivative of the mean-attractive energy with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cMeanAttractiveEnergyCSW_d2Volume        ! Second derivative of the mean-attractive energy with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cEffectiveRDF_d2Volume                   ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cEffectiveRDFYukawa_d2Volume             ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cEffectiveRDFSutherland_d2Volume         ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cEffectiveRDFCSW_d2Volume                ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPackingFraction_d2Volume              ! Second derivative of the effective packing fraction with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPackingFractionYK_d2Volume            ! Second derivative of the effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPackingFractionSL_d2Volume            ! Second derivative of the effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPackingFractionCSW_d2Volume           ! Second derivative of the effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPFraction_d2Volume_dDensity           ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPFractionYK_d2Volume_dDensity         ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPFractionSL_d2Volume_dDensity         ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffPFractionCSW_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cRDFunction_d2Volume                     ! Second derivative of the contact radial distribution function with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2cMeanEnergyFluctuations_d2Volume         ! Second derivative of the mean-attractive energy fluctuations with respect to the volume
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2FirstTPTCoeff_d2Volume_dDensity          ! First derivative with respect to the volume of the cross derivative of the first-order perturbation coefficient with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffectiveRDF_d2Volume_dDensity           ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffectiveRDFYK_d2Volume_dDensity         ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffectiveRDFCSW_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: d2EffectiveRDFSL_d2Volume_dDensity         ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: d2AuxEffRDF_d2Volume_dDensity              ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: d2AuxEffRDFYukawa_d2Volume_dDensity        ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: d2AuxEffRDFSutherland_d2Volume_dDensity    ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: d2AuxEffRDFCSW_d2Volume_dDensity           ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
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
      IF( PYHCBCorrectionLogical ) CSW2CoefficientsMatrix = CSWAllCoefficientA2
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
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
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
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
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
!            FIRST DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE VOLUME             !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! First derivative of the number density with respect to the volume [mol / (Å³ . m³)]
dNumberDensity_dVolume = - rNumberDensity / mVolume

! First derivative of the inverse number density with respect to the volume [mol . (Å³ / m³)]
dInverseNumberDensity_dVolume = 1.D0 / rNumberDensity / mVolume

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dIdealFreeEnergy_dVolume = - (1.D0 / mVolume) * cUniversalGas * Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! First derivative of the reduced densities with respect to the volume
drDensity_dVolume(0) = rDensityConstants(0) * dNumberDensity_dVolume ! [mol / (Å³ . m³)]
drDensity_dVolume(1) = rDensityConstants(1) * dNumberDensity_dVolume ! [mol / (Å² . m³)]
drDensity_dVolume(2) = rDensityConstants(2) * dNumberDensity_dVolume ! [mol / (Å . m³)]
drDensity_dVolume(3) = rDensityConstants(3) * dNumberDensity_dVolume ! [mol / m³]

! Auxiliary factors of the derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume ! [mol / (Å³ . m³)]
dAuxiliaryHSFactor_dVolume(1) = - ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * &
&     ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * drDensity_dVolume(3) + DLOG( 1.D0 - rDensity(3) ) * ( ( 3.D0 * rDensity(2) * rDensity(2) &
&     * drDensity_dVolume(2) / ( rDensity(3) * rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * &
&     drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) * rDensity(3) ) ) - drDensity_dVolume(0) )
dAuxiliaryHSFactor_dVolume(2) = 3.D0 * ( drDensity_dVolume(1) * rDensity(2) + drDensity_dVolume(2) * rDensity(1) ) / ( 1.D0 - &
&     rDensity(3) ) + 3.D0 * rDensity(1) * rDensity(2) * drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) )
dAuxiliaryHSFactor_dVolume(3) = ( rDensity(2) * rDensity(2) * ( 3.D0 * drDensity_dVolume(2) * rDensity(3) * rDensity(3) - &
&     ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * drDensity_dVolume(2) ) * rDensity(3) + drDensity_dVolume(3) * &
&     rDensity(2) ) ) / ( ( rDensity(3) - 1.D0 ) * ( rDensity(3) - 1.D0 ) * ( rDensity(3) - 1.D0 ) * rDensity(3) * rDensity(3) )

! First derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume [mol / m³]
dHSBoublikMansoori_dVolume = HSBoublikMansoori * rNumberDensity * dInverseNumberDensity_dVolume + ( 6.D0 / (cPi * &
&     rNumberDensity) ) * SUM( dAuxiliaryHSFactor_dVolume )

! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dExcludedVolumeFreeEnergy_dVolume = dHSBoublikMansoori_dVolume * AngleAverageMixture * cUniversalGas * Temperature ! Proven units

! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume (Boublik's term) [Joule / m³ or Pa]
IF( ReferenceBoublikLogical ) THEN
  dHCBBoublik_dVolume = - ( ( ( sumSqRadius * sumArea * sumArea ) / ( 9.D0 * sumVolume * sumVolume ) ) - 1.D0 ) * ( 1.D0 / &
  &     ( 1.D0 - rDensity(3) ) ) + ( ( sumRadius * sumArea ) / ( sumVolume * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + &
  &     ( ( sumSqRadius * sumArea * sumArea * ( 1.D0 + rDensity(3) ) ) / ( 9.D0 * sumVolume * sumVolume * ( 1.D0 - rDensity(3) ) * &
  &     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
  dHCBBoublik_dVolume = dHCBBoublik_dVolume * drDensity_dVolume(3)
  ! Excluded-volume contribution to the Helmholtz free energy (molar basis) [J / mol]
  dExcludedVolumeFreeEnergy_dVolume = dHCBBoublik_dVolume * cUniversalGas * Temperature ! Proven units
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! First derivative of the mean-attractive energy with respect to the volume
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! First derivative of the effective packing fraction with respect to the volume [mol / m³]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) &
        &     + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) &
        &     + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        dEffPackingFractionSL_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensity(3) + 2.D0 * EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) &
        &     + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        dEffPackingFractionYK_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensity(3) + 2.D0 * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) &
        &     + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent)
        dEffPackingFractionCSW_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * &
        &     rDensity(3) + 2.D0 * EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + 3.D0 * &
        &     rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsCSW(3,iComponent,jComponent)
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture &
        &     + 3.D0 * rDensityMixture * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        dEffPackingFractionSL_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        dEffPackingFractionYK_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dEffPackingFraction_dVolume(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture &
        &     + 3.D0 * rDensityMixture * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent)
        dEffPackingFractionCSW_dVolume(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture + 3.D0 * rDensityMixture * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent)
      END IF
    END IF
    dEffPackingFraction_dVolume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume
    IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      dEffPackingFractionSL_dVolume(iComponent,jComponent) = - dEffPackingFractionSL_dVolume(iComponent,jComponent) / mVolume
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      dEffPackingFractionYK_dVolume(iComponent,jComponent) = - dEffPackingFractionYK_dVolume(iComponent,jComponent) / mVolume
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
      dEffPackingFractionCSW_dVolume(iComponent,jComponent) = - dEffPackingFractionCSW_dVolume(iComponent,jComponent) / mVolume
    END IF
    ! Cross derivative of the effective packing fraction with respect to the volume and the density [mol . (Å³ / m³)]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * &
        &     rDensity(3) ) + ( 9.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent) * &
        &     rDensity(3) ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * &
        &     rDensity(3) ) )
        dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) = dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / &
        &     mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( &
        &     EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * rDensity(3) ) + ( 4.D0 * &
        &     EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensity(3) * rDensity(3) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * &
        &     rDensity(3) ) )
        dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) = dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / &
        &     mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( &
        &     EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * rDensity(3) ) + ( 4.D0 * &
        &     EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensity(3) * rDensity(3) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = - ( ( 2.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) ) + ( 6.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) ) ) / ( mVolume * rNumberDensity )
        dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) = - ( ( 2.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) ) + ( 6.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) ) ) / ( mVolume * rNumberDensity )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture ) + ( 9.D0 * rDensityMixture * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture ) )
        dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) = dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / &
        &     mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( &
        &     EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 4.D0 * EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / mVolume - &
        &     ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture ) )
        dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) = dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / &
        &     mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( &
        &     EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * rDensityMixture ) + ( 4.D0 * &
        &     EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensityMixture * rDensityMixture ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dxEffPFraction_dVolume_dDensity(iComponent,jComponent) = - ( ( 2.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) ) + ( 6.D0 * rDensityMixture * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent) ) ) / ( mVolume * rNumberDensity )
        dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) = - ( ( 2.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) ) + ( 6.D0 * rDensityMixture * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) ) ) / ( mVolume * rNumberDensity )
      END IF
    END IF
    ! First derivative of the contact radial distribution function with respect to the volume [mol / m³]
    dcRDFunction_dVolume(iComponent,jComponent) = ( drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) * &
    &     ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
    &     rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) + ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * &
    &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * &
    &     rDensity(3) / ( 1.D0 - rDensity(3) ) ) + ( 3.D0 * rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
    &     rDensity(3) ) ) ) ) ) )
    ! First derivative of the radial distribution function for an effective packing fraction with respect to the volume [mol / m³]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dcEffectiveRDF_dVolume(iComponent,jComponent) = ( dEffPackingFraction_dVolume(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) * ( 1.D0 + &
        &     ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
        &     EffPackingFraction(iComponent,jComponent) / ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( 2.D0 * &
        &     AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * EffPackingFraction(iComponent,jComponent) / ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) ) + ( 3.D0 * EffPackingFraction(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) ) ) )
      ELSE
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dcEffectiveRDF_dVolume(iComponent,jComponent) = dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) = ( dEffPackingFractionSL_dVolume(iComponent,jComponent) / ( ( &
        &     1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     ) ) + ( 2.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     ) + ( 3.D0 * EffPackingFractionSutherland(iComponent,jComponent) * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) = ( dEffPackingFractionYK_dVolume(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) * &
        &     ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) / ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * EffPackingFractionYukawa(iComponent,jComponent) / ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 3.D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) ) &
        &     ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dcEffectiveRDFCSW_dVolume(iComponent,jComponent) = dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dcEffectiveRDF_dVolume(iComponent,jComponent) = - ( 0.5D0 * dEffPackingFraction_dVolume(iComponent,jComponent) / ( ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( EffPackingFraction(iComponent,jComponent) / &
        &     2.D0 ) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) * dEffPackingFraction_dVolume(iComponent,jComponent)
      ELSE
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        dcEffectiveRDF_dVolume(iComponent,jComponent) = ( dEffPackingFraction_dVolume(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + &
        &     ( 3.D0 * gAux1 * ( dEffPackingFraction_dVolume(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 2.D0 &
        &     * gAux2 * EffPackingFraction(iComponent,jComponent) * ( ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) = - ( 0.5D0 * dEffPackingFractionSL_dVolume(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) + ( ( 3.D0 * ( 1.D0 - ( EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) * dEffPackingFractionSL_dVolume(iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) = - ( 0.5D0 * dEffPackingFractionYK_dVolume(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( ( 3.D0 * ( 1.D0 - ( EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) ) * dEffPackingFractionYK_dVolume(iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        dcEffectiveRDFCSW_dVolume(iComponent,jComponent) = ( dEffPackingFractionCSW_dVolume(iComponent,jComponent) / &
        &     ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) ) ) + ( 3.D0 * gAux1 * ( dEffPackingFractionCSW_dVolume(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) &
        &     ) + ( ( 2.D0 * gAux2 * EffPackingFractionCSW(iComponent,jComponent) * ( ( 2.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * dEffPackingFractionCSW_dVolume(iComponent,jComponent) ) ) / ( ( &
        &     1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) ) )
      END IF
    END IF
    ! First derivative of the mean-attractive energy with respect to the volume [K . mol / m³]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) = - 2.D0 * cPi * ( ijPotentialRangeCubic(iComponent,jComponent) - &
      &     1.D0 ) * ( ijaDiameterSphereCubic(iComponent,jComponent) / 3.D0 ) * ( ( dNumberDensity_dVolume * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + dcEffectiveRDF_dVolume(iComponent,jComponent) * &
      &     rNumberDensity ) * ijaWellDepth(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) * &
      &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + &
      &     ( rNumberDensity * dcEffectiveRDF_dVolume(iComponent,jComponent) ) )
      dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) * &
      &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) + &
      &     ( rNumberDensity * dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) ) ) * &
      &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + &
      &     ( rNumberDensity * dcEffectiveRDF_dVolume(iComponent,jComponent) ) )
      dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) = - cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) * &
      &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) + &
      &     ( rNumberDensity * dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) = - ijaWellDepth(iComponent,jComponent) * ( &
      &     ijSecondVirialCoefficientField(iComponent,jComponent) - ijSecondVirialCoefficient(iComponent,jComponent) ) * ( &
      &     dNumberDensity_dVolume * cEffectiveRadialDistributionFunction(iComponent,jComponent) + rNumberDensity * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) )
      dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) = - ijaWellDepth(iComponent,jComponent) * ( &
      &     ijSecondVirialCoefficientField(iComponent,jComponent) - ijSecondVirialCoefficient(iComponent,jComponent) ) * ( &
      &     dNumberDensity_dVolume * cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) + rNumberDensity * &
      &     dcEffectiveRDFCSW_dVolume(iComponent,jComponent) )
    END IF
  END DO
END DO

! Mixing rule for the derivative of the first-order perturbation coefficient with respect to the volume [K . mol / m³]
CALL Mixing_Rules( mFraction, dcMeanAttractiveEnergy_dVolume, dpFirstOrderCoefficient_dVolume )

! First derivative of the first-order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMeanAttractiveFreeEnergy_dVolume = dpFirstOrderCoefficient_dVolume * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the hard-sphere isothermal compressibility with respect to the volume [mol / m³]
dHSIsothermalCompressibility_dVolume = ( ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
&     * drDensity_dVolume(3) * rDensity(0) ) ) / ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * &
&     rDensity(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ) ) - &
&     ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) / &
&     ( ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ) * ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * &
&     rDensity(2) ) ) ) ) * ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * rDensity(0) * drDensity_dVolume(3) ) + ( 6.D0 * ( ( drDensity_dVolume(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( drDensity_dVolume(2) * rDensity(1) * ( 1.D0 - rDensity(3) ) ) - ( drDensity_dVolume(3) * rDensity(1) * &
&     rDensity(2) ) ) ) + ( 27.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) )

! First derivative of the hard-convex-body isothermal compressibility with respect to the volume [mol / m³]
IF( PYHCBCorrectionLogical ) THEN
  HCBAux1 = sumRadius * sumArea / sumVolume - 1.D0
  HCBAux2 = 1.D0 - 2.D0 * sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume )
  dHCBIsothermalCompressibility_dVolume = - ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
  &     drDensity_dVolume(3) / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) * ( 1.D0 + 2.D0 * &
  &     rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + &
  &     2.D0 * rDensity(3) * HCBAux2 + 2.D0 * rDensity(3) * rDensity(3) * HCBAux2 )
ELSE
  HCBAux1 = sumRadius * sumArea / sumVolume - 1.D0
  HCBAux2 = 1.D0 - 2.D0 * sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume )
  HCBAux3 = sumSqRadius * sumArea * sumArea / ( 9.D0 * sumVolume * sumVolume )
  dHCBIsothermalCompressibility_dVolume = - ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
  &     ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * rDensity(3) * &
  &     rDensity(3) * HCBAux2 - 12.D0 * rDensity(3) * rDensity(3) * HCBAux3 ) * drDensity_dVolume(3) / ( ( 1.D0 + 2.D0 * &
  &     rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     rDensity(3) - 4.D0 * HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + &
  &     rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * HCBAux3 * &
  &     rDensity(3) * rDensity(3) * rDensity(3) ) )
END IF

! First derivative of the Zhang's factor with respect to the number of particles [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    dZhangFactor_dVolume = 4.D0 * ZhangCorrection * ( cPi * SUM( mFraction * aDiameterSphere * aDiameterSphere * aDiameterSphere ) &
    &     / 6.D0 ) * ( cPi * SUM( mFraction * aDiameterSphere * aDiameterSphere * aDiameterSphere ) / 6.D0 ) * rNumberDensity * &
    &     dNumberDensity_dVolume
  ELSE
    dZhangFactor_dVolume = - 0.25D0 * ZhangCorrection * SecondVirialCoefficientMixtureSingle * &
    &     SecondVirialCoefficientMixtureSingle * rNumberDensity * rNumberDensity / mVolume
  END IF
END IF

! First derivative of the mean-attractive energy fluctuation with respect to the volume
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Cross derivative of the effective radial distribution function with respect to the volume and the density [mol . (Å³ / m³)]
    IF( EffPFractionMixingRule == 1 ) THEN
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) = ( ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) + ( 2.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + &
        &     ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * &
        &     ( ( ( ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) &
        &     ) ) + ( dEffPackingFraction_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 3.D0 * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) ) ) + ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( &
        &     2.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( &
        &     1.D0 + EffPackingFraction(iComponent,jComponent) ) ) + ( ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     EffPackingFraction(iComponent,jComponent) * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 4.D0 * &
        &     ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * EffPackingFraction(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) ) )
      ELSE
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 &
        &     + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * &
        &     EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) + 4.D0 * ( 1.D0 + dAux1 + 2.D0 * &
        &     EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + EffPackingFraction(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) = ( ( &
        &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) + ( 2.D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     ) ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) + ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * ( ( ( ( &
        &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 3.D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) ) + ( 2.D0 * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) &
        &     * ( 1.D0 + EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( ( 2.D0 + &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * EffPackingFractionSutherland(iComponent,jComponent) &
        &     * dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 4.D0 * ( 2.D0 + &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * EffPackingFractionSutherland(iComponent,jComponent) * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) &
        &     ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) = ( ( &
        &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) + ( 2.D0 * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( 3.D0 * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * ( ( ( ( &
        &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) + ( dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 3.D0 * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) ) ) ) + ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( &
        &     2.D0 * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) + ( ( 2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) * EffPackingFractionYukawa(iComponent,jComponent) &
        &     * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) ) + ( ( 4.D0 * ( 2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     EffPackingFractionYukawa(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) = dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * &
        &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) ) + dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dDensity(iComponent,jComponent) * ( ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) &
        &     + 4.D0 * ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) - ( ( 3.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 3.D0 * &
        &     ( 1.D0 - ( EffPackingFraction(iComponent,jComponent) / 2.D0 ) ) * &
        &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 12.D0 * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) / 2.D0 ) * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) )
      ELSE
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )         
        dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) = ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + &
        &     ( ( dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( &
        &     3.D0 * gAux1 + 2.D0 ) + 3.D0 * gAux1 * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + &
        &     ( ( ( dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     ( ( 9.D0 * gAux1 + 4.D0 * gAux2 ) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) ) ) + 2.D0 * gAux2 * &
        &     EffPackingFraction(iComponent,jComponent) * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * ( 2.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 8.D0 * gAux2 * EffPackingFraction(iComponent,jComponent) * ( &
        &     2.D0 + EffPackingFraction(iComponent,jComponent) ) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) ) * &
        &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 12.D0 * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) ) * &
        &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) ) + ( ( 12.D0 * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )         
        dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) = ( dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) &
        &     / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) ) ) + ( ( dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dDensity(iComponent,jComponent) * ( 3.D0 * gAux1 + 2.D0 ) + 3.D0 * gAux1 * &
        &     dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) &
        &     ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) ) + ( ( &
        &     ( dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dDensity(iComponent,jComponent) * &
        &     ( ( 9.D0 * gAux1 + 4.D0 * gAux2 ) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) ) ) + 2.D0 * gAux2 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * ( 2.D0 &
        &     + EffPackingFractionCSW(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) ) + ( ( 8.D0 * gAux2 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dDensity(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) )
      END IF
    END IF
    ! Cross derivative of the first-order perturbation coefficient with respect to the volume and the density [K . mol . (Å³ / m³)]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) = ( dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( cMeanAttractiveEnergy(iComponent,jComponent) * dNumberDensity_dVolume / &
      &     ( rNumberDensity * rNumberDensity ) ) + ( ( ( dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunct_dDensity(iComponent,jComponent) ) - ( ( cMeanAttractiveEnergy(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) * dcEffectiveRDF_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunct_dDensity(iComponent,jComponent) ) + ( cMeanAttractiveEnergy(iComponent,jComponent) * &
      &     dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) ) ) * ( 1.D0 / &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) = ( &
      &     dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) / rNumberDensity ) - &
      &     ( cMeanAttractiveEnergySutherland(iComponent,jComponent) * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity &
      &     ) ) + ( ( ( dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) ) - ( ( &
      &     cMeanAttractiveEnergySutherland(iComponent,jComponent) / cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) &
      &     ) * dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) ) + ( &
      &     cMeanAttractiveEnergySutherland(iComponent,jComponent) * &
      &     dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) ) ) * ( 1.D0 / &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) = ( &
      &     dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) / rNumberDensity ) - &
      &     ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) &
      &     + ( ( ( dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) ) - ( ( &
      &     cMeanAttractiveEnergyYukawa(iComponent,jComponent) / cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) * &
      &     dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) * dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) ) &
      &     + ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) &
      &     ) ) * ( 1.D0 / cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) = ( dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * dNumberDensity_dVolume / ( rNumberDensity * &
      &     rNumberDensity ) ) + ( ( ( dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) ) - ( ( cMeanAttractiveEnergyCSW(iComponent,jComponent) &
      &     / cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) * dcEffectiveRDFCSW_dVolume(iComponent,jComponent) &
      &     * dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) ) + ( cMeanAttractiveEnergyCSW(iComponent,jComponent) &
      &     * dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) ) ) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) )
    END IF
    ! First derivative of the mean-attractive energy fluctuation with respect to the volume [K² . mol / m³]
    IF( ZhangCorrectionLogical ) THEN
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * rNumberDensity &
        &     * ZhangFactor ) + ( dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * HSIsothermalCompressibility * &
        &     ZhangFactor ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) &
        &     * rNumberDensity * ZhangFactor ) + ( rNumberDensity * HSIsothermalCompressibility * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dZhangFactor_dVolume ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * rNumberDensity &
        &     * ZhangFactor ) + ( dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
        &     HSIsothermalCompressibility * ZhangFactor ) + ( HSIsothermalCompressibility * &
        &     dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity * ZhangFactor ) + &
        &     ( rNumberDensity * HSIsothermalCompressibility * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
        &     dZhangFactor_dVolume ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * rNumberDensity &
        &     * ZhangFactor ) + ( dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * HSIsothermalCompressibility &
        &     * ZhangFactor ) + ( HSIsothermalCompressibility * &
        &     dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity * ZhangFactor ) + &
        &     ( rNumberDensity * HSIsothermalCompressibility * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
        &     dZhangFactor_dVolume ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHCBIsothermalCompressibility_dVolume * rNumberDensity * &
        &     ZhangFactor ) + ( dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * HCBIsothermalCompressibility * &
        &     ZhangFactor ) + ( HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * ZhangFactor ) + ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * rNumberDensity * &
        &     HCBIsothermalCompressibility * dZhangFactor_dVolume ) )
      END IF
    ELSE
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * &
        &     rNumberDensity ) + ( dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * HSIsothermalCompressibility ) + &
        &     ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * &
        &     rNumberDensity ) + ( dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * HSIsothermalCompressibility ) &
        &     + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * &
        &     ( ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHSIsothermalCompressibility_dVolume * &
        &     rNumberDensity ) + ( dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * HSIsothermalCompressibility ) + &
        &     ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
        dcMeanEnergyFluctuations_dVolume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dHCBIsothermalCompressibility_dVolume * rNumberDensity ) + ( &
        &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * HCBIsothermalCompressibility ) + ( &
        &     HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity ) )
      END IF
    END IF
  END DO
END DO

! Mixing rule for the second-order perturbation coefficient [K² . mol / m³]
CALL Mixing_Rules( mFraction, dcMeanEnergyFluctuations_dVolume, dpSecondOrderCoefficient_dVolume )

! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMeanFluctuationFreeEnergy_dVolume = dpSecondOrderCoefficient_dVolume * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  dHigherOrderTerms_dVolume = 0.D0
  DO nOrder = 3, nHigherOrder
    dHigherOrderTerms_dVolume = dHigherOrderTerms_dVolume + ( dMeanAttractiveFreeEnergy_dVolume * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) ) + ( ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy &
    &     * dMeanAttractiveFreeEnergy_dVolume) ) / ( mMeanAttractiveHelmholtzFreeEnergy ) ) )
  END DO
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
  dHigherOrderFEnergy_dVolume = dHigherOrderTerms_dVolume ! Proven units
ELSE
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
  dHigherOrderFEnergy_dVolume = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! First derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMonomerFreeEnergy_dVolume = dExcludedVolumeFreeEnergy_dVolume + dMeanAttractiveFreeEnergy_dVolume + &
&     dMeanFluctuationFreeEnergy_dVolume + dHigherOrderFEnergy_dVolume ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! First derivative of the total Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dTotalFreeEnergy_dVolume = dIdealFreeEnergy_dVolume + dMonomerFreeEnergy_dVolume ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!            SECOND DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE VOLUME            !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Second derivative of the number density with respect to the volume [mol² / (Å³ . m⁶)]
d2NumberDensity_d2Volume = 2.D0 * rNumberDensity / ( mVolume * mVolume )

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2IdealFreeEnergy_d2Volume = ( 1.D0 / ( mVolume * mVolume ) ) * cUniversalGas * Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Second derivative of the reduced densities with respect to the volume
d2rDensity_d2Volume(0) = rDensityConstants(0) * d2NumberDensity_d2Volume ! [mol² / (Å³ . m⁶)]
d2rDensity_d2Volume(1) = rDensityConstants(1) * d2NumberDensity_d2Volume ! [mol² / (Å² . m⁶)]
d2rDensity_d2Volume(2) = rDensityConstants(2) * d2NumberDensity_d2Volume ! [mol² / (Å . m⁶)]
d2rDensity_d2Volume(3) = rDensityConstants(3) * d2NumberDensity_d2Volume ! [mol² / m⁶]

! Auxiliary factors of the derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume ! [mol² / (Å³ . m⁶)]
d2AuxiliaryHSFactor_d2Volume(1) = - ( ( ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) &
&     * rDensity(3) ) ) - drDensity_dVolume(0) ) * ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * drDensity_dVolume(3) ) - ( ( ( rDensity(2) &
&     * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * drDensity_dVolume(3) * drDensity_dVolume(3) &
&     / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * &
&     rDensity(3) ) ) - rDensity(0) ) * ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * d2rDensity_d2Volume(3) ) - ( ( 1.D0 / ( 1.D0 - &
&     rDensity(3) ) ) * drDensity_dVolume(3) * ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) &
&     * rDensity(3) ) ) - drDensity_dVolume(0) ) ) + ( DLOG( 1.D0 - rDensity(3) ) * ( ( ( 3.D0 * ( ( 2.D0 * rDensity(2) * &
&     drDensity_dVolume(2) * drDensity_dVolume(2) ) + ( rDensity(2) * rDensity(2) * d2rDensity_d2Volume(2) ) ) ) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 6.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) * drDensity_dVolume(3) ) / ( rDensity(3) * &
&     rDensity(3) * rDensity(3) ) ) - ( ( 2.D0 * ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) * &
&     drDensity_dVolume(3) ) + ( d2rDensity_d2Volume(3) * rDensity(2) * rDensity(2) * rDensity(2) ) ) ) / ( rDensity(3) * &
&     rDensity(3) * rDensity(3) ) ) + ( ( 6.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) ) ) - d2rDensity_d2Volume(0) ) )
d2AuxiliaryHSFactor_d2Volume(2) = ( ( 3.D0 * ( ( d2rDensity_d2Volume(1) * rDensity(2) ) + ( 2.D0 * drDensity_dVolume(1) * &
&     drDensity_dVolume(2) ) + ( d2rDensity_d2Volume(2) * rDensity(1) ) ) ) / ( 1.D0 - rDensity(3) ) ) + ( ( 3.D0 * &
&     ( ( drDensity_dVolume(1) * rDensity(2) ) + ( drDensity_dVolume(2) * rDensity(1) ) ) * drDensity_dVolume(3) ) / ( ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 3.D0 * ( ( drDensity_dVolume(1) * drDensity_dVolume(3) * rDensity(2) ) + &
&     ( drDensity_dVolume(2) * drDensity_dVolume(3) * rDensity(1) ) + ( rDensity(1) * rDensity(2) * d2rDensity_d2Volume(3) ) ) ) / &
&     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 6.D0 * rDensity(1) * rDensity(2) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
d2AuxiliaryHSFactor_d2Volume(3) = - ( ( ( 2.D0 * rDensity(2) * drDensity_dVolume(2) * ( 3.D0 * drDensity_dVolume(2) * rDensity(3) &
&     * rDensity(3) - ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * drDensity_dVolume(2) ) * rDensity(3) + &
&     drDensity_dVolume(3) * rDensity(2) ) ) + ( rDensity(2) * rDensity(2) * ( ( 6.D0 * rDensity(3) * drDensity_dVolume(2) * &
&     drDensity_dVolume(3) ) + ( 3.D0 * rDensity(3) * rDensity(3) * d2rDensity_d2Volume(2) ) + ( d2rDensity_d2Volume(3) * &
&     rDensity(2) ) + ( drDensity_dVolume(3) * drDensity_dVolume(2) ) - ( drDensity_dVolume(3) * ( 3.D0 * drDensity_dVolume(3) * &
&     rDensity(2) + 3.D0 * drDensity_dVolume(2) ) ) - ( 3.D0 * rDensity(3) * ( ( d2rDensity_d2Volume(3) * rDensity(2) ) + &
&     ( drDensity_dVolume(3) * drDensity_dVolume(2) ) + d2rDensity_d2Volume(2) ) ) ) ) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) * rDensity(3) ) ) + ( ( ( rDensity(2) * rDensity(2) * ( 3.D0 * &
&     drDensity_dVolume(2) * rDensity(3) * rDensity(3) - ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * &
&     drDensity_dVolume(2) ) * rDensity(3) + drDensity_dVolume(3) * rDensity(2) ) ) * ( drDensity_dVolume(3) * ( 2.D0 - ( 5.D0 * &
&     rDensity(3) ) ) ) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     rDensity(3) * rDensity(3) * rDensity(3) ) )

! Second derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume [mol² / m⁶]
d2HSBoublikMansoori_d2Volume = ( 12.D0 / cPi ) * dInverseNumberDensity_dVolume * SUM( dAuxiliaryHSFactor_dVolume ) + &
&     ( 6.D0 / ( cPi * rNumberDensity ) ) * SUM( d2AuxiliaryHSFactor_d2Volume )

! Second derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2ExcludedVolumeFreeEnergy_d2Volume = d2HSBoublikMansoori_d2Volume * AngleAverageMixture * cUniversalGas * Temperature ! Proven units

! Second derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume (Boublik's term) [Pa . mol / m³]
IF( ReferenceBoublikLogical ) THEN
  d2HCBBoublik_d2Volume = - ( ( ( sumSqRadius * sumArea * sumArea ) / ( 9.D0 * sumVolume * sumVolume ) ) - 1.D0 ) * ( 1.D0 / &
  &     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 2.D0 * sumRadius * sumArea ) / ( sumVolume * ( 1.D0 - &
  &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( sumSqRadius * sumArea * sumArea * ( 4.D0 + &
  &     2.D0 * rDensity(3) ) ) / ( 9.D0 * sumVolume * sumVolume * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
  d2HCBBoublik_d2Volume = d2HCBBoublik_d2Volume * drDensity_dVolume(3) * drDensity_dVolume(3)
  d2HCBBoublik_d2Volume = d2HCBBoublik_d2Volume + dHCBBoublik_dVolume * d2rDensity_d2Volume(3) / drDensity_dVolume(3)
  ! Excluded-volume contribution to the Helmholtz free energy (molar basis) [J / mol]
  d2ExcludedVolumeFreeEnergy_d2Volume = d2HCBBoublik_d2Volume * cUniversalGas * Temperature ! Proven units
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! Second derivative of the mean-attractive energy with respect to the volume
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Second derivative of the effective packing fraction with respect to the volume [mol² / m⁶]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) + 9.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) ) / ( mVolume * mVolume )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) ) / ( mVolume * mVolume )
        d2EffPackingFractionSL_d2Volume(iComponent,jComponent) = - dEffPackingFractionSL_dVolume(iComponent,jComponent) / mVolume &
        &     + ( EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * rDensity(3) + 4.D0 * rDensity(3) * &
        &     rDensity(3) * EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) ) / ( mVolume * mVolume )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) ) / ( mVolume * mVolume )
        d2EffPackingFractionYK_d2Volume(iComponent,jComponent) = - dEffPackingFractionYK_dVolume(iComponent,jComponent) / mVolume &
        &     + ( EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) ) / ( mVolume * mVolume )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = ( ( 2.D0 * EffPackingFractionCoefficients(1,iComponent,jComponent) &
        &     * rDensity(3) ) + ( 6.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) ) + ( &
        &     12.D0 * EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( &
        &     mVolume * mVolume )
        d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) = ( ( 2.D0 * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) ) + ( 6.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) ) + ( 12.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) ) ) / ( mVolume * mVolume )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + 4.D0 * &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture + 9.D0 * rDensityMixture &
        &     * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent) ) / ( mVolume * &
        &     mVolume )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + 4.D0 * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficients(2,iComponent,jComponent) ) / ( mVolume * mVolume )
        d2EffPackingFractionSL_d2Volume(iComponent,jComponent) = - dEffPackingFractionSL_dVolume(iComponent,jComponent) / mVolume &
        &     + ( EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * rDensityMixture + 4.D0 * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) ) / ( mVolume * mVolume )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = - dEffPackingFraction_dVolume(iComponent,jComponent) / mVolume + &
        &     ( EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + 4.D0 * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficients(2,iComponent,jComponent) ) / ( mVolume * mVolume )
        d2EffPackingFractionYK_d2Volume(iComponent,jComponent) = - dEffPackingFractionYK_dVolume(iComponent,jComponent) / mVolume &
        &     + ( EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * rDensityMixture + 4.D0 * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) ) / ( mVolume * mVolume )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2EffPackingFraction_d2Volume(iComponent,jComponent) = ( ( 2.D0 * EffPackingFractionCoefficients(1,iComponent,jComponent) &
        &     * rDensityMixture ) + ( 6.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture ) + ( 12.D0 * EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture * rDensityMixture ) ) / ( mVolume * mVolume )
        d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) = ( ( 2.D0 * rDensityMixture * &
        &     EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) ) + ( 6.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) ) + ( 12.D0 * rDensityMixture * rDensityMixture * &
        &     rDensityMixture * EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) ) ) / ( mVolume * mVolume )
      END IF
    END IF
    ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density [mol² . (Å³ / m⁶)]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) ) + ( 8.D0 * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2,iComponent,jComponent) ) + ( 27.D0 * rDensity(3) * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) ) + ( 8.D0 * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2,iComponent,jComponent) ) )
        d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) = dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) / &
        &     mVolume - dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( &
        &     rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensity(3) ) + ( 8.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) ) + ( 8.D0 * &
        &     rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2,iComponent,jComponent) ) )
        d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) = dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) / &
        &     mVolume - dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( &
        &     rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensity(3) ) + ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) &
        &     * rDensity(3) * rDensity(3) / ( rNumberDensity * mVolume * mVolume ) ) + ( 18.D0 * rDensity(3) * rDensity(3) * &
        &     rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent) / ( rNumberDensity * mVolume * mVolume ) )
        d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) = ( 4.D0 * rDensity(3) * rDensity(3) * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) / ( rNumberDensity * mVolume * mVolume ) ) + ( 18.D0 * &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3) / ( &
        &     rNumberDensity * mVolume * mVolume ) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture ) + ( 8.D0 * &
        &     rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(2,iComponent,jComponent) ) + ( 27.D0 * &
        &     rDensityMixture * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture ) + ( 8.D0 * &
        &     rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(2,iComponent,jComponent) ) )
        d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) = dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) / &
        &     mVolume - dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( &
        &     rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 8.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = dxEffPFraction_dVolume_dDensity(iComponent,jComponent) / mVolume &
        &     - dEffPackingFraction_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * &
        &     mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture ) + ( 8.D0 * &
        &     rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(2,iComponent,jComponent) ) )
        d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) = dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) / &
        &     mVolume - dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / ( mVolume * mVolume ) + ( 1.D0 / ( &
        &     rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensityMixture ) + ( 8.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) = ( 4.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) &
        &     * rDensityMixture * rDensityMixture / ( rNumberDensity * mVolume * mVolume ) ) + ( 18.D0 * rDensityMixture * &
        &     rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent) / ( rNumberDensity * &
        &     mVolume * mVolume ) )
        d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) = ( 4.D0 * rDensityMixture * rDensityMixture * &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) / ( rNumberDensity * mVolume * mVolume ) ) + ( 18.D0 * &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensityMixture * rDensityMixture * rDensityMixture / ( &
        &     rNumberDensity * mVolume * mVolume ) )
      END IF
    END IF
    ! Second derivative of the contact radial distribution function with respect to the volume [mol² / m⁶]
    d2cRDFunction_d2Volume(iComponent,jComponent) = ( ( d2rDensity_d2Volume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
    &     rDensity(3) ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
    &     ( 1.D0 + 2.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) + ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * &
    &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * &
    &     rDensity(3) / ( 1.D0 - rDensity(3) ) ) + ( 3.D0 * rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
    &     rDensity(3) ) ) ) ) ) ) ) + ( ( 2.D0 * drDensity_dVolume(3) * drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
    &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) * ( 1.D0 + ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
    &     AuxDiameterRelationship * ( 2.D0 + ( 3.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) ) + ( 2.D0 * &
    &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
    &     AuxDiameterRelationship * AuxDiameterRelationship * ( 1.D0 + ( 6.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * &
    &     rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) ) ) ) )
    ! Second derivative of the radial distribution function for an effective packing fraction with respect to the volume [mol² / m⁶]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        d2cEffectiveRDF_d2Volume(iComponent,jComponent) = ( ( d2EffPackingFraction_d2Volume(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) * ( 1.D0 + &
        &     ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
        &     EffPackingFraction(iComponent,jComponent) / ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( 2.D0 * &
        &     AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * EffPackingFraction(iComponent,jComponent) / ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) ) + ( 3.D0 * EffPackingFraction(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) ) ) ) ) + ( ( ( 2.D0 * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) * ( 1.D0 + ( ( 3.D0 * AuxDiameterRelationship * &
        &     cAverageDiameterRelationship(iComponent,jComponent) ) * ( 2.D0 + ( 3.D0 * EffPackingFraction(iComponent,jComponent) &
        &     / ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) ) + ( ( 2.D0 * AuxDiameterRelationship * &
        &     AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) ) * ( 1.D0 + ( 6.D0 * EffPackingFraction(iComponent,jComponent) &
        &     / ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( 6.D0 * EffPackingFraction(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) ) ) ) )
      ELSE
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        d2cEffectiveRDF_d2Volume(iComponent,jComponent) = d2EffPackingFraction_d2Volume(iComponent,jComponent) * ( 1.D0 + dAux1 &
        &     + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + EffPackingFraction(iComponent,jComponent) &
        &     * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * ( ( 2.D0 &
        &     * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + 4.D0 * ( &
        &     1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2cEffectiveRDFSutherland_d2Volume(iComponent,jComponent) = ( ( d2EffPackingFractionSL_d2Volume(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + 2.D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     ) ) + ( 2.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     ) + ( 3.D0 * EffPackingFractionSutherland(iComponent,jComponent) * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) ) ) ) + ( ( ( 2.D0 * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) * ( 1.D0 + ( ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) ) * ( 2.D0 &
        &     + ( 3.D0 * EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) + ( ( 2.D0 * AuxDiameterRelationship * &
        &     AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) ) * ( 1.D0 + ( 6.D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     ) + ( 6.D0 * EffPackingFractionSutherland(iComponent,jComponent) * &
        &     EffPackingFractionSutherland(iComponent,jComponent) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2cEffectiveRDFYukawa_d2Volume(iComponent,jComponent) = ( ( d2EffPackingFractionYK_d2Volume(iComponent,jComponent) / ( ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( 1.D0 + &
        &     2.D0 * EffPackingFractionYukawa(iComponent,jComponent) / ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) &
        &     ) ) + ( 2.D0 * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * ( ( 2.D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) / ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( &
        &     3.D0 * EffPackingFractionYukawa(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) / ( ( 1.D0 &
        &     - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) &
        &     ) ) ) ) + ( ( ( 2.D0 * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) ) ) * ( 1.D0 + ( ( 3.D0 * AuxDiameterRelationship * &
        &     cAverageDiameterRelationship(iComponent,jComponent) ) * ( 2.D0 + ( 3.D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) / ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) ) + &
        &     ( ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) ) * ( 1.D0 + ( 6.D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) / ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( &
        &     6.D0 * EffPackingFractionYukawa(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) / ( ( 1.D0 &
        &     - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) &
        &     ) ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        d2cEffectiveRDFCSW_d2Volume(iComponent,jComponent) = d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) * ( 1.D0 + &
        &     dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) ) + dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * ( ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) + &
        &     4.D0 * ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        d2cEffectiveRDF_d2Volume(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     d2EffPackingFraction_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) - ( ( 3.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 3.D0 * &
        &     ( 1.D0 - ( EffPackingFraction(iComponent,jComponent) / 2.D0 ) ) * &
        &     d2EffPackingFraction_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 12.D0 * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) / 2.D0 ) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) )
      ELSE
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        d2cEffectiveRDF_d2Volume(iComponent,jComponent) = d2EffPackingFraction_d2Volume(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( ( &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * ( 2.D0 + &
        &     3.D0 * gAux1 ) ) + ( 3.D0 * gAux1 * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     d2EffPackingFraction_d2Volume(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) &
        &     * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + &
        &     ( ( ( 4.D0 * gAux2 + 9.D0 * gAux1 ) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) + &
        &     ( 2.D0 * gAux2 * EffPackingFraction(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) &
        &     * d2EffPackingFraction_d2Volume(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( 8.D0 * gAux2 * &
        &     EffPackingFraction(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2cEffectiveRDFSutherland_d2Volume(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) ) * &
        &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 12.D0 * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2cEffectiveRDFYukawa_d2Volume(iComponent,jComponent) = - ( ( 0.5D0 * &
        &     d2EffPackingFractionYK_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) ) * d2EffPackingFractionYK_d2Volume(iComponent,jComponent) &
        &     ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 12.D0 * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        d2cEffectiveRDFCSW_d2Volume(iComponent,jComponent) = d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) / ( ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) + ( ( &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * ( &
        &     2.D0 + 3.D0 * gAux1 ) ) + ( 3.D0 * gAux1 * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) ) + ( ( ( 4.D0 * gAux2 + 9.D0 * gAux1 ) * ( 1.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) ) + ( 2.D0 * gAux2 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) + ( &
        &     8.D0 * gAux2 * EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) )
      END IF
    END IF
    ! Second derivative of the mean-attractive energy with respect to the volume [K . mol² / m⁶]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) = - 2.D0 * cPi * ( ijPotentialRangeCubic(iComponent,jComponent) - &
      &     1.D0 ) * ( ijaDiameterSphereCubic(iComponent,jComponent) / 3.D0 ) * ( ( d2NumberDensity_d2Volume * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + ( 2.D0 * dcEffectiveRDF_dVolume(iComponent,jComponent) &
      &     * dNumberDensity_dVolume ) + ( rNumberDensity * d2cEffectiveRDF_d2Volume(iComponent,jComponent) ) ) * &
      &     ijaWellDepth(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) * ( ( &
      &     d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) * dNumberDensity_dVolume ) + ( rNumberDensity * &
      &     d2cEffectiveRDF_d2Volume(iComponent,jComponent) ) )
      d2cMeanAttractiveEnergySL_d2Volume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) * &
      &     ( ( d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) * dNumberDensity_dVolume ) + ( rNumberDensity * &
      &     d2cEffectiveRDFSutherland_d2Volume(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) = - 2.D0 * cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) ) ) * ( ( &
      &     d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunction(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) * dNumberDensity_dVolume ) + ( rNumberDensity * &
      &     d2cEffectiveRDF_d2Volume(iComponent,jComponent) ) )
      d2cMeanAttractiveEnergyYK_d2Volume(iComponent,jComponent) = - cPi * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) * ( ( &
      &     d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) * dNumberDensity_dVolume ) + ( rNumberDensity * &
      &     d2cEffectiveRDFYukawa_d2Volume(iComponent,jComponent) ) )
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) = - ijaWellDepth(iComponent,jComponent) * ( &
      &     ijSecondVirialCoefficientField(iComponent,jComponent) - ijSecondVirialCoefficient(iComponent,jComponent) ) * ( ( &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) * d2NumberDensity_d2Volume ) + ( 2.D0 * &
      &     dNumberDensity_dVolume * dcEffectiveRDF_dVolume(iComponent,jComponent) ) + ( &
      &     d2cEffectiveRDF_d2Volume(iComponent,jComponent) * rNumberDensity ) )
      d2cMeanAttractiveEnergyCSW_d2Volume(iComponent,jComponent) = - ijaWellDepth(iComponent,jComponent) * ( &
      &     ijSecondVirialCoefficientField(iComponent,jComponent) - ijSecondVirialCoefficient(iComponent,jComponent) ) * ( ( &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) * d2NumberDensity_d2Volume ) + ( 2.D0 * &
      &     dNumberDensity_dVolume * dcEffectiveRDFCSW_dVolume(iComponent,jComponent) ) + ( rNumberDensity * &
      &     d2cEffectiveRDFCSW_d2Volume(iComponent,jComponent) ) )
    END IF
  END DO
END DO

! Mixing rule for the derivative of the first-order perturbation coefficient with respect to the volume [K . mol² / m⁶]
CALL Mixing_Rules( mFraction, d2cMeanAttractiveEnergy_d2Volume, d2pFirstOrderCoefficient_d2Volume )

! Second derivative of the first-order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MeanAttractiveFreeEnergy_d2Volume = d2pFirstOrderCoefficient_d2Volume * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Auxiliary factors of a second derivative of the hard-sphere isothermal compressibility with respect to the volume
d2AuxHSIC_d2Volume(1) = ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) ) - ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     drDensity_dVolume(3) * rDensity(0) ) ) ! [mol / (Å³ . m³)]
d2AuxHSIC_d2Volume(2) = ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * &
&     ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ! [1 / Å³]
d2AuxHSIC_d2Volume(3) = ( d2rDensity_d2Volume(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) ) ) - ( 8.D0 * drDensity_dVolume(0) * drDensity_dVolume(3) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( 12.D0 * rDensity(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) * &
&     drDensity_dVolume(3) * drDensity_dVolume(3) ) - ( 4.D0 * rDensity(0) * d2rDensity_d2Volume(3) * ( ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) ! [mol² / (Å³ . m⁶)]
d2AuxHSIC_d2Volume(4) = ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * rDensity(0) * drDensity_dVolume(3) ) + ( 6.D0 * ( ( drDensity_dVolume(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( drDensity_dVolume(2) * rDensity(1) * ( 1.D0 - rDensity(3) ) ) - ( drDensity_dVolume(3) * rDensity(1) * &
&     rDensity(2) ) ) ) + ( 27.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) ! [mol / (Å³ . m³)]
d2AuxHSIC_d2Volume(5) = ( d2AuxHSIC_d2Volume(1) / ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) - ( ( ( 2.D0 * rDensity(0) * &
&     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) / &
&     ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * d2AuxHSIC_d2Volume(4) ) ! [mol . (Å³ / m³)]
d2AuxHSIC_d2Volume(6) = - ( ( rDensity(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) ) / ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * ( ( d2rDensity_d2Volume(0) * ( ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( 4.D0 * drDensity_dVolume(0) * drDensity_dVolume(3) * ( 1.D0 - rDensity(3) ) ) &
&     - ( 2.D0 * rDensity(0) * ( 1.D0 - rDensity(3) ) * d2rDensity_d2Volume(3) ) + ( 2.D0 * rDensity(0) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) + ( 6.D0 * ( ( d2rDensity_d2Volume(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * drDensity_dVolume(1) * drDensity_dVolume(2) ) - ( 2.D0 * drDensity_dVolume(1) * drDensity_dVolume(3) * &
&     rDensity(2) ) + ( rDensity(1) * d2rDensity_d2Volume(2) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * rDensity(1) * &
&     drDensity_dVolume(2) * drDensity_dVolume(3) ) - ( rDensity(1) * rDensity(2) * d2rDensity_d2Volume(3) ) ) ) + ( 54.D0 * &
&     rDensity(2) * drDensity_dVolume(2) * drDensity_dVolume(2) ) + ( 27.D0 * rDensity(2) * rDensity(2) * &
&     d2rDensity_d2Volume(2) ) ) ! [mol² / m⁶]

! Second derivative of the hard-sphere isothermal compressibility with respect to the volume [mol² / m⁶]
d2HSIsothermalCompressibility_d2Volume = ( d2AuxHSIC_d2Volume(3) / d2AuxHSIC_d2Volume(2) ) - ( ( d2AuxHSIC_d2Volume(1) / &
&     ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * d2AuxHSIC_d2Volume(4) ) - ( d2AuxHSIC_d2Volume(5) * &
&     d2AuxHSIC_d2Volume(4) ) + d2AuxHSIC_d2Volume(6)

! Second derivative of the hard-convex-body isothermal compressibility with respect to the volume [mol² / m⁶]
IF( PYHCBCorrectionLogical ) THEN
  HCBAux1 = sumRadius * sumArea / sumVolume - 1.D0
  HCBAux2 = 1.D0 - 2.D0 * sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume )
  d2HCBIsothermalCompressibility_d2Volume = - ( d2rDensity_d2Volume(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 &
  &     - rDensity(3) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * &
  &     rDensity(3) * rDensity(3) * HCBAux2 ) ) / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) &
  &     * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) ) - drDensity_dVolume(3) * &
  &     drDensity_dVolume(3) * ( ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 6.D0 * HCBAux1 + &
  &     2.D0 * HCBAux2 + 4.D0 * rDensity(3) * HCBAux2 ) - ( 3.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 4.D0 + 2.D0 &
  &     * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * rDensity(3) * rDensity(3) * HCBAux2 ) ) ) &
  &     / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 &
  &     + rDensity(3) * rDensity(3) * HCBAux2 ) ) - ( 2.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * rDensity(3) &
  &     * rDensity(3) * HCBAux2 ) * ( 2.D0 * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 ) ) / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 &
  &     + rDensity(3) * rDensity(3) * HCBAux2 ) * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) * &
  &     ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 ) ) )
ELSE
  HCBAux1 = sumRadius * sumArea / sumVolume - 1.D0
  HCBAux2 = 1.D0 - 2.D0 * sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume )
  HCBAux3 = sumSqRadius * sumArea * sumArea / ( 9.D0 * sumVolume * sumVolume )
  d2HCBIsothermalCompressibility_d2Volume = - ( d2rDensity_d2Volume(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 &
  &     - rDensity(3) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * &
  &     rDensity(3) * rDensity(3) * HCBAux2 - 12.D0 * rDensity(3) * rDensity(3) * HCBAux3 ) ) / ( ( 1.D0 + 2.D0 * rDensity(3) * &
  &     HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * &
  &     HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * &
  &     HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * HCBAux3 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) ) ) - drDensity_dVolume(3) * drDensity_dVolume(3) * ( ( ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( &
  &     1.D0 - rDensity(3) ) * ( 6.D0 * HCBAux1 + 2.D0 * HCBAux2 + 4.D0 * rDensity(3) * HCBAux2 - 24.D0 * rDensity(3) * HCBAux3 ) &
  &     ) - ( 3.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + &
  &     2.D0 * rDensity(3) * HCBAux2 + 2.D0 * rDensity(3) * rDensity(3) * HCBAux2 - 12.D0 * rDensity(3) * rDensity(3) * HCBAux3 ) &
  &     ) ) / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) &
  &     * rDensity(3) * rDensity(3) - 4.D0 * HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) * ( 1.D0 + 2.D0 * rDensity(3) * &
  &     HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * &
  &     HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) ) - ( 2.D0 * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( &
  &     1.D0 - rDensity(3) ) ) * ( 4.D0 + 2.D0 * HCBAux1 + 6.D0 * rDensity(3) * HCBAux1 + 2.D0 * rDensity(3) * HCBAux2 + 2.D0 * &
  &     rDensity(3) * rDensity(3) * HCBAux2 - 12.D0 * rDensity(3) * rDensity(3) * HCBAux3 ) * ( 2.D0 * HCBAux1 + 2.D0 * &
  &     rDensity(3) * HCBAux2 + 4.D0 * rDensity(3) * rDensity(3) * rDensity(3) * HCBAux3 - 12.D0 * rDensity(3) * rDensity(3) * &
  &     HCBAux3 ) ) / ( ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * &
  &     rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) * ( 1.D0 + 2.D0 * &
  &     rDensity(3) * HCBAux1 + rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     rDensity(3) - 4.D0 * HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) ) * ( 1.D0 + 2.D0 * rDensity(3) * HCBAux1 + &
  &     rDensity(3) * rDensity(3) * HCBAux2 + HCBAux3 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * HCBAux3 * &
  &     rDensity(3) * rDensity(3) * rDensity(3) ) ) )
END IF

! Second derivative of the Zhang's factor with respect to the number of particles [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    d2ZhangFactor_d2Volume = 4.D0 * ZhangCorrection * ( cPi * SUM( mFraction * aDiameterSphere * aDiameterSphere * &
    &     aDiameterSphere ) / 6.D0 ) * ( cPi * SUM( mFraction * aDiameterSphere * aDiameterSphere * aDiameterSphere ) / 6.D0 ) * &
    &     ( ( dNumberDensity_dVolume * dNumberDensity_dVolume ) + ( rNumberDensity * d2NumberDensity_d2Volume ) )
  ELSE
    d2ZhangFactor_d2Volume = 0.75D0 * ZhangCorrection * SecondVirialCoefficientMixtureSingle * &
    &     SecondVirialCoefficientMixtureSingle * rNumberDensity * rNumberDensity / ( mVolume * mVolume )
  END IF
END IF

! Second derivative of the mean-attractive energy fluctuation with respect to the volume
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! First derivative with respect to the volume of an auxiliary factor of the cross derivative of the effective radial distribution function with respect to the volume and the density [mol² . (Å³ / m⁶)]
    IF( .NOT. PotentialTypeLogical(4) ) THEN
      d2AuxEffRDF_d2Volume_dDensity(1,iComponent,jComponent) = ( ( d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) * ( &
      &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) ) + ( 6.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) + ( 2.D0 * &
      &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( ( dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) + ( dEffPackingFraction_dDensity(iComponent,jComponent) * &
      &     d2EffPackingFraction_d2Volume(iComponent,jComponent) ) ) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 &
      &     - EffPackingFraction(iComponent,jComponent) ) )
      d2AuxEffRDF_d2Volume_dDensity(2,iComponent,jComponent) = 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
      &     AuxDiameterRelationship * ( ( ( ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) * &
      &     EffPackingFraction(iComponent,jComponent) ) * d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) ) - ( 2.D0 * &
      &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) ) + ( 2.D0 * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) &
      &     * dEffPackingFraction_dVolume(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) ) + ( 2.D0 &
      &     * d2EffPackingFraction_d2Volume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( 2.D0 &
      &     + EffPackingFraction(iComponent,jComponent) ) ) + ( 2.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 &
      &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( ( 4.D0 * &
      &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) ) ) + ( 8.D0 * &
      &     dEffPackingFraction_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) ) ) / ( ( &
      &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) ) ) )
      d2AuxEffRDF_d2Volume_dDensity(3,iComponent,jComponent) = 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
      &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 &
      &     * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * ( 1.D0 &
      &     + 4.D0 * EffPackingFraction(iComponent,jComponent) + EffPackingFraction(iComponent,jComponent) * &
      &     EffPackingFraction(iComponent,jComponent) ) ) + ( 2.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * &
      &     d2EffPackingFraction_d2Volume(iComponent,jComponent) * ( 1.D0 + 4.D0 * EffPackingFraction(iComponent,jComponent) + &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) ) ) + ( 4.D0 * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     dEffPackingFraction_dDensity(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) ) + &
      &     ( d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 2.D0 - &
      &     EffPackingFraction(iComponent,jComponent) - EffPackingFraction(iComponent,jComponent) * &
      &     EffPackingFraction(iComponent,jComponent) ) ) + ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * ( 2.D0 - EffPackingFraction(iComponent,jComponent) - &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) ) ) - &
      &     ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * ( 1.D0 + 2.D0 * EffPackingFraction(iComponent,jComponent) ) ) ) / &
      &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 &
      &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 5.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     ( ( 2.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * &
      &     ( 1.D0 + 4.D0 * EffPackingFraction(iComponent,jComponent) + EffPackingFraction(iComponent,jComponent) * &
      &     EffPackingFraction(iComponent,jComponent) ) ) + ( dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * &
      &     EffPackingFraction(iComponent,jComponent) * ( 2.D0 - EffPackingFraction(iComponent,jComponent) - &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) ) ) ) ) / ( ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) )
    ELSE
      geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
      &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
      dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
      &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
      &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
      dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
      &     sumVolume * geoAux )
      d2AuxEffRDF_d2Volume_dDensity(1,iComponent,jComponent) = d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) * ( 1.D0 &
      &     + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / &
      &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + &
      &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     ( ( ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( &
      &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( 4.D0 * ( &
      &     1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / &
      &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) ) ) )
      d2AuxEffRDF_d2Volume_dDensity(2,iComponent,jComponent) = ( d2EffPackingFraction_d2Volume(iComponent,jComponent) * &
      &     dEffPackingFraction_dDensity(iComponent,jComponent) + dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) * ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
      &     ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) + 4.D0 * &
      &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) &
      &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
      &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
      &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) )
      d2AuxEffRDF_d2Volume_dDensity(3,iComponent,jComponent) = dEffPackingFraction_dVolume(iComponent,jComponent) * &
      &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( ( ( &
      &     3.D0 * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) + &
      &     2.D0 * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( 5.D0 * ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
      &     ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) + 4.D0 * &
      &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) &
      &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
      &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 &
      &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) )
    END IF
    IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      d2AuxEffRDFSutherland_d2Volume_dDensity(1,iComponent,jComponent) = ( ( &
      &     d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( 2.D0 * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) ) + ( 6.D0 * &
      &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) ) + ( 2.D0 * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( ( dEffPackingFractionSL_dVolume(iComponent,jComponent) * &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) + ( &
      &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * d2EffPackingFractionSL_d2Volume(iComponent,jComponent) &
      &     ) ) ) ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) )
      d2AuxEffRDFSutherland_d2Volume_dDensity(2,iComponent,jComponent) = 3.D0 * &
      &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * ( ( ( ( ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) ) * &
      &     d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) ) - ( 2.D0 * &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) ) + ( 2.D0 * &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( &
      &     2.D0 + EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) &
      &     * ( 2.D0 + EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( ( 4.D0 * &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( &
      &     1.D0 - EffPackingFractionSutherland(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) ) ) + &
      &     ( 8.D0 * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( 2.D0 + &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     ) ) )
      d2AuxEffRDFSutherland_d2Volume_dDensity(3,iComponent,jComponent) = 2.D0 * &
      &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
      &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 * &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( &
      &     1.D0 + 4.D0 * EffPackingFractionSutherland(iComponent,jComponent) + &
      &     EffPackingFractionSutherland(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( 2.D0 &
      &     * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
      &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) * ( 1.D0 + 4.D0 * &
      &     EffPackingFractionSutherland(iComponent,jComponent) + EffPackingFractionSutherland(iComponent,jComponent) * &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( 4.D0 * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * ( 2.D0 + &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) &
      &     * EffPackingFractionSutherland(iComponent,jComponent) * ( 2.D0 - EffPackingFractionSutherland(iComponent,jComponent) - &
      &     EffPackingFractionSutherland(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( &
      &     2.D0 - EffPackingFractionSutherland(iComponent,jComponent) - EffPackingFractionSutherland(iComponent,jComponent) * &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) - ( dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * &
      &     EffPackingFractionSutherland(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( 1.D0 + &
      &     2.D0 * EffPackingFractionSutherland(iComponent,jComponent) ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     ) ) + ( ( 5.D0 * dEffPackingFractionSL_dVolume(iComponent,jComponent) * ( ( 2.D0 * &
      &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
      &     ( 1.D0 + 4.D0 * EffPackingFractionSutherland(iComponent,jComponent) + &
      &     EffPackingFractionSutherland(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( &
      &     dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) * ( &
      &     2.D0 - EffPackingFractionSutherland(iComponent,jComponent) - EffPackingFractionSutherland(iComponent,jComponent) * &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) )
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      d2AuxEffRDFYukawa_d2Volume_dDensity(1,iComponent,jComponent) = ( ( d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) &
      &     * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 2.D0 * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) + ( 6.D0 * &
      &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) + ( 2.D0 * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( ( dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
      &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) + ( &
      &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * d2EffPackingFractionYK_d2Volume(iComponent,jComponent) ) ) &
      &     ) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) )
      d2AuxEffRDFYukawa_d2Volume_dDensity(2,iComponent,jComponent) = 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
      &     AuxDiameterRelationship * ( ( ( ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) ) - ( &
      &     2.D0 * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) + ( 2.D0 * &
      &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( &
      &     2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     d2EffPackingFractionYK_d2Volume(iComponent,jComponent) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * ( &
      &     2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
      &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
      &     ( ( ( 4.D0 * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 8.D0 * &
      &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( 2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) ) ) &
      &     / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) &
      &     ) ) )
      d2AuxEffRDFYukawa_d2Volume_dDensity(3,iComponent,jComponent) = 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
      &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 &
      &     * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( &
      &     1.D0 + 4.D0 * EffPackingFractionYukawa(iComponent,jComponent) + EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 2.D0 * &
      &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * d2EffPackingFractionYK_d2Volume(iComponent,jComponent) * &
      &     ( 1.D0 + 4.D0 * EffPackingFractionYukawa(iComponent,jComponent) + EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( 4.D0 * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * ( &
      &     2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( &
      &     d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) * ( 2.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) - EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( 2.D0 - EffPackingFractionYukawa(iComponent,jComponent) - &
      &     EffPackingFractionYukawa(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) ) ) - &
      &     ( dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( 1.D0 + 2.D0 * EffPackingFractionYukawa(iComponent,jComponent) &
      &     ) ) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) &
      &     ) ) + ( ( 5.D0 * dEffPackingFractionYK_dVolume(iComponent,jComponent) * ( ( 2.D0 * &
      &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * &
      &     ( 1.D0 + 4.D0 * EffPackingFractionYukawa(iComponent,jComponent) + EffPackingFractionYukawa(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) * &
      &     EffPackingFractionYukawa(iComponent,jComponent) * ( 2.D0 - EffPackingFractionYukawa(iComponent,jComponent) - &
      &     EffPackingFractionYukawa(iComponent,jComponent) * EffPackingFractionYukawa(iComponent,jComponent) ) ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
      &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) &
      &     * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionYukawa(iComponent,jComponent) ) ) ) )
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
      &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
      dAux1 = ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
      &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
      &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( sumVolume * geoAux )
      dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
      &     sumVolume * geoAux )
      d2AuxEffRDFCSW_d2Volume_dDensity(1,iComponent,jComponent) = d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) * ( &
      &     1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 &
      &     ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) ) + dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * &
      &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * ( ( ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) ) + ( &
      &     4.D0 * ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
      &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 &
      &     ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) ) )
      d2AuxEffRDFCSW_d2Volume_dDensity(2,iComponent,jComponent) = ( d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) * &
      &     dEffPackingFractionCSW_dDensity(iComponent,jComponent) + dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
      &     dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) ) * ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) + 4.D0 * ( 1.D0 + dAux1 + 2.D0 * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + EffPackingFractionCSW(iComponent,jComponent) * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) )
      d2AuxEffRDFCSW_d2Volume_dDensity(3,iComponent,jComponent) = dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
      &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dDensity(iComponent,jComponent) * ( ( &
      &     ( 3.D0 * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) &
      &     ) + 2.D0 * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - dAux1 + dAux2 ) ) / ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) ) ) + ( 5.D0 * ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 2.D0 * ( dAux2 - 1.D0 ) + 2.D0 * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) + 4.D0 * ( 1.D0 + dAux1 + 2.D0 * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + EffPackingFractionCSW(iComponent,jComponent) * &
      &     EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) ) / ( ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
      &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) ) )
    END IF
    ! First derivative with respect to the volume of the cross derivative of the effective radial distribution function with respect to the volume and the density [mol² . (Å³ / m⁶)]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      d2EffectiveRDF_d2Volume_dDensity(iComponent,jComponent) = SUM( d2AuxEffRDF_d2Volume_dDensity(:,iComponent,jComponent) )
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffectiveRDFSL_d2Volume_dDensity(iComponent,jComponent) = SUM( &
        &     d2AuxEffRDFSutherland_d2Volume_dDensity(:,iComponent,jComponent) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffectiveRDFYK_d2Volume_dDensity(iComponent,jComponent) = SUM( &
        &     d2AuxEffRDFYukawa_d2Volume_dDensity(:,iComponent,jComponent) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2EffectiveRDFCSW_d2Volume_dDensity(iComponent,jComponent) = SUM( &
        &     d2AuxEffRDFCSW_d2Volume_dDensity(:,iComponent,jComponent) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        d2EffectiveRDF_d2Volume_dDensity(iComponent,jComponent) = -0.5D0 * d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( ( 3.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFraction(iComponent,jComponent) ) * d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) * d2EffPackingFraction_d2Volume(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) - ( ( &
        &     6.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( &
        &     12.D0 * ( 1.D0 - 0.5D0 * EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) * d2EffPackingFraction_d2Volume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 24.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFraction(iComponent,jComponent) ) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dxEffPFraction_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) - ( ( &
        &     18.D0 * dEffPackingFraction_dDensity(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 60.D0 * &
        &     ( 1.D0 - 0.5D0 * EffPackingFraction(iComponent,jComponent) ) * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dVolume(iComponent,jComponent) ) / ( ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) )
      ELSE
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        d2AuxEffRDF_d2Volume_dDensity(1,iComponent,jComponent) = d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) + ( ( &
        &     3.D0 * gAux1 + 2.D0 ) * ( 2.D0 * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) + dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFraction_d2Volume(iComponent,jComponent) ) + 3.D0 * gAux1 * &
        &     d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) ) / &
        &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) )
        d2AuxEffRDF_d2Volume_dDensity(2,iComponent,jComponent) = ( 2.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( 9.D0 * &
        &     gAux1 + 2.D0 * gAux2 + 3.D0 ) + dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFraction_d2Volume(iComponent,jComponent) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) * ( 9.D0 &
        &     * gAux1 + 4.D0 * gAux2 ) + 2.D0 * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) * ( 9.D0 * &
        &     gAux1 + 4.D0 * gAux2 ) + 2.D0 * d2EffPFraction_d2Volume_dDensity(iComponent,jComponent) * gAux2 * &
        &     EffPackingFraction(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) )
        d2AuxEffRDF_d2Volume_dDensity(3,iComponent,jComponent) = ( 16.D0 * dxEffPFraction_dVolume_dDensity(iComponent,jComponent) &
        &     * dEffPackingFraction_dVolume(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 2.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) * gAux2 + 4.D0 * dEffPackingFraction_dVolume(iComponent,jComponent) * &
        &     dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 9.D0 * gAux1 + 8.D0 * gAux2 ) + 8.D0 * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) * d2EffPackingFraction_d2Volume(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * gAux2 ) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) + ( 40.D0 * gAux2 * dEffPackingFraction_dVolume(iComponent,jComponent) &
        &     * dEffPackingFraction_dVolume(iComponent,jComponent) * dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) * ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) )
        d2EffectiveRDF_d2Volume_dDensity(iComponent,jComponent) = SUM( d2AuxEffRDF_d2Volume_dDensity(:,iComponent,jComponent) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2EffectiveRDFSL_d2Volume_dDensity(iComponent,jComponent) = - 0.5D0 * &
        &     d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( ( 3.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * d2EffPFractionSL_d2Volume_dDensity(iComponent,jComponent) ) &
        &     / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) - ( ( 6.D0 * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 12.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionSL_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) + ( ( 24.D0 * ( 1.D0 - 0.5D0 * EffPackingFractionSutherland(iComponent,jComponent) ) * &
        &     dEffPackingFractionSL_dVolume(iComponent,jComponent) * dxEffPFractionSL_dVolume_dDensity(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) - ( ( 18.D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) + ( ( 60.D0 * ( 1.D0 - 0.5D0 * EffPackingFractionSutherland(iComponent,jComponent) ) * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * dEffPackingFractionSL_dVolume(iComponent,jComponent) &
        &     * dEffPackingFractionSL_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2EffectiveRDFYK_d2Volume_dDensity(iComponent,jComponent) = - 0.5D0 * &
        &     d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( ( 3.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * d2EffPFractionYK_d2Volume_dDensity(iComponent,jComponent) ) / ( &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) - ( ( 3.D0 * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * d2EffPackingFractionYK_d2Volume(iComponent,jComponent) ) &
        &     / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) - ( ( 6.D0 * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) / &
        &     ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 12.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionYK_d2Volume(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 24.D0 * ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dxEffPFractionYK_dVolume_dDensity(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) - &
        &     ( ( 18.D0 * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( ( 60.D0 * ( 1.D0 - 0.5D0 * EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * dEffPackingFractionYK_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionYK_dVolume(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        d2AuxEffRDFCSW_d2Volume_dDensity(1,iComponent,jComponent) = d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) &
        &     + ( ( 3.D0 * gAux1 + 2.D0 ) * ( 2.D0 * dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) + dEffPackingFractionCSW_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) ) + 3.D0 * gAux1 * &
        &     d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) &
        &     ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
        d2AuxEffRDFCSW_d2Volume_dDensity(2,iComponent,jComponent) = ( 2.D0 * dEffPackingFractionCSW_dVolume(iComponent,jComponent) &
        &     * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dDensity(iComponent,jComponent) * ( &
        &     9.D0 * gAux1 + 2.D0 * gAux2 + 3.D0 ) + dEffPackingFractionCSW_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 9.D0 * gAux1 + 4.D0 * gAux2 ) + 2.D0 * dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * ( &
        &     9.D0 * gAux1 + 4.D0 * gAux2 ) + 2.D0 * d2EffPFractionCSW_d2Volume_dDensity(iComponent,jComponent) * gAux2 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + EffPackingFractionCSW(iComponent,jComponent) ) ) / ( ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
        d2AuxEffRDFCSW_d2Volume_dDensity(3,iComponent,jComponent) = ( 16.D0 * EffPackingFractionCSW(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dxEffPFractionCSW_dVolume_dDensity(iComponent,jComponent) * &
        &     ( 2.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * gAux2 + 4.D0 * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionCSW(iComponent,jComponent) ) * ( &
        &     9.D0 * gAux1 + 8.D0 * gAux2 ) + 8.D0 * dEffPackingFractionCSW_dDensity(iComponent,jComponent) * &
        &     d2EffPackingFractionCSW_d2Volume(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * gAux2 ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) &
        &     * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) &
        &     * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) &
        &     ) + ( 40.D0 * gAux2 * dEffPackingFractionCSW_dVolume(iComponent,jComponent) * &
        &     dEffPackingFractionCSW_dVolume(iComponent,jComponent) * dEffPackingFractionCSW_dDensity(iComponent,jComponent) * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 2.D0 + EffPackingFractionCSW(iComponent,jComponent) ) ) / ( ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
        d2EffectiveRDFCSW_d2Volume_dDensity(iComponent,jComponent) = SUM( &
        &     d2AuxEffRDFCSW_d2Volume_dDensity(:,iComponent,jComponent) )
      END IF
    END IF
    ! First derivative with respect to the volume of the cross derivative of the first-order perturbation coefficient with respect to the volume and the density [K . mol² . (Å³ / m⁶)]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) = ( d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( 2.D0 * dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * dNumberDensity_dVolume / &
      &     ( rNumberDensity * rNumberDensity ) ) - ( cMeanAttractiveEnergy(iComponent,jComponent) * d2NumberDensity_d2Volume / &
      &     ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * cMeanAttractiveEnergy(iComponent,jComponent) * dNumberDensity_dVolume &
      &     * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( 1.D0 / &
      &     ( cEffectiveRadialDistributionFunction(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) ) * dcEffectiveRDF_dVolume(iComponent,jComponent) * ( ( &
      &     dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * dEffRadialDistributionFunct_dDensity(iComponent,jComponent) &
      &     ) - ( ( cMeanAttractiveEnergy(iComponent,jComponent) / cEffectiveRadialDistributionFunction(iComponent,jComponent) ) * &
      &     dEffRadialDistributionFunct_dDensity(iComponent,jComponent) * dcEffectiveRDF_dVolume(iComponent,jComponent) ) + &
      &     ( cMeanAttractiveEnergy(iComponent,jComponent) * dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) ) ) ) + ( ( &
      &     1.D0 / cEffectiveRadialDistributionFunction(iComponent,jComponent) ) * ( ( &
      &     d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) * dEffRadialDistributionFunct_dDensity(iComponent,jComponent) &
      &     ) + ( 2.D0 * dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * &
      &     dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) ) - ( ( 1.D0 / &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) * ( ( &
      &     dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * dEffRadialDistributionFunct_dDensity(iComponent,jComponent) * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) ) + ( cMeanAttractiveEnergy(iComponent,jComponent) * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) * dxEffectiveRDF_dVolume_dDensity(iComponent,jComponent) ) + &
      &     ( dEffRadialDistributionFunct_dDensity(iComponent,jComponent) * cMeanAttractiveEnergy(iComponent,jComponent) * &
      &     d2cEffectiveRDF_d2Volume(iComponent,jComponent) ) ) ) + ( ( 1.D0 / ( &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) ) ) * ( cMeanAttractiveEnergy(iComponent,jComponent) * &
      &     dEffRadialDistributionFunct_dDensity(iComponent,jComponent) * dcEffectiveRDF_dVolume(iComponent,jComponent) * &
      &     dcEffectiveRDF_dVolume(iComponent,jComponent) ) ) + ( cMeanAttractiveEnergy(iComponent,jComponent) * &
      &     d2EffectiveRDF_d2Volume_dDensity(iComponent,jComponent) ) ) )
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) = ( d2cMeanAttractiveEnergySL_d2Volume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( 2.D0 * dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * dNumberDensity_dVolume / &
      &     ( rNumberDensity * rNumberDensity ) ) - ( cMeanAttractiveEnergySutherland(iComponent,jComponent) * &
      &     d2NumberDensity_d2Volume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * &
      &     cMeanAttractiveEnergySutherland(iComponent,jComponent) * dNumberDensity_dVolume * dNumberDensity_dVolume / ( &
      &     rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( 1.D0 / ( &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) ) * &
      &     dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) * ( ( &
      &     dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) ) - ( ( &
      &     cMeanAttractiveEnergySutherland(iComponent,jComponent) / cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) &
      &     ) * dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) * &
      &     dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) ) + ( cMeanAttractiveEnergySutherland(iComponent,jComponent) * &
      &     dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) ) ) ) + ( ( 1.D0 / &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) * ( ( &
      &     d2cMeanAttractiveEnergySL_d2Volume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
      &     dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) ) - ( ( 1.D0 / &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) * ( ( &
      &     dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) * &
      &     dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) ) + ( cMeanAttractiveEnergySutherland(iComponent,jComponent) * &
      &     dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) * &
      &     dxEffectiveRDFSutherland_dVolume_dDensity(iComponent,jComponent) ) + ( &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) * cMeanAttractiveEnergySutherland(iComponent,jComponent) &
      &     * d2cEffectiveRDFSutherland_d2Volume(iComponent,jComponent) ) ) ) + ( ( 1.D0 / ( &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) ) ) * ( &
      &     cMeanAttractiveEnergySutherland(iComponent,jComponent) * dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) &
      &     * dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) * dcEffectiveRDFSutherland_dVolume(iComponent,jComponent) ) &
      &     ) + ( cMeanAttractiveEnergySutherland(iComponent,jComponent) * &
      &     d2EffectiveRDFSL_d2Volume_dDensity(iComponent,jComponent) ) ) )
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) = ( d2cMeanAttractiveEnergyYK_d2Volume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( 2.D0 * dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * dNumberDensity_dVolume / &
      &     ( rNumberDensity * rNumberDensity ) ) - ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * &
      &     d2NumberDensity_d2Volume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * &
      &     cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dNumberDensity_dVolume * dNumberDensity_dVolume / ( &
      &     rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( 1.D0 / ( &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) ) * &
      &     dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) * ( ( dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) &
      &     * dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) ) - ( ( &
      &     cMeanAttractiveEnergyYukawa(iComponent,jComponent) / cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) * &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) * dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) ) &
      &     + ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) &
      &     ) ) ) + ( ( 1.D0 / cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) * ( ( &
      &     d2cMeanAttractiveEnergyYK_d2Volume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) ) + ( 2.D0 * &
      &     dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * &
      &     dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) ) - ( ( 1.D0 / &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) * ( ( &
      &     dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) * dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) ) &
      &     + ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) * &
      &     dxEffectiveRDFYukawa_dVolume_dDensity(iComponent,jComponent) ) + ( &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) * cMeanAttractiveEnergyYukawa(iComponent,jComponent) * &
      &     d2cEffectiveRDFYukawa_d2Volume(iComponent,jComponent) ) ) ) + ( ( 1.D0 / ( &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) ) ) * ( &
      &     cMeanAttractiveEnergyYukawa(iComponent,jComponent) * dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) * &
      &     dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) * dcEffectiveRDFYukawa_dVolume(iComponent,jComponent) ) ) + &
      &     ( cMeanAttractiveEnergyYukawa(iComponent,jComponent) * d2EffectiveRDFYK_d2Volume_dDensity(iComponent,jComponent) ) ) )
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) = ( d2cMeanAttractiveEnergyCSW_d2Volume(iComponent,jComponent) / &
      &     rNumberDensity ) - ( dNumberDensity_dVolume * dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) / ( &
      &     rNumberDensity * rNumberDensity ) ) + ( 2.D0 * dNumberDensity_dVolume * dNumberDensity_dVolume * &
      &     cMeanAttractiveEnergyCSW(iComponent,jComponent) / ( rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( &
      &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * dNumberDensity_dVolume + &
      &     cMeanAttractiveEnergyCSW(iComponent,jComponent) * d2NumberDensity_d2Volume ) / ( rNumberDensity * rNumberDensity ) ) + &
      &     ( ( ( d2cMeanAttractiveEnergyCSW_d2Volume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) ) + ( &
      &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) ) &
      &     - ( ( dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * dcEffectiveRDFCSW_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) + ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * &
      &     d2cEffectiveRDFCSW_d2Volume(iComponent,jComponent) * dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) + ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * &
      &     dcEffectiveRDFCSW_dVolume(iComponent,jComponent) * dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) - ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * &
      &     dcEffectiveRDFCSW_dVolume(iComponent,jComponent) * dcEffectiveRDFCSW_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) / ( &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) ) ) + ( &
      &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) ) &
      &     + ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * d2EffectiveRDFCSW_d2Volume_dDensity(iComponent,jComponent) ) ) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) - ( ( ( ( &
      &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) ) - ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * &
      &     dcEffectiveRDFCSW_dVolume(iComponent,jComponent) * dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) + ( cMeanAttractiveEnergyCSW(iComponent,jComponent) * &
      &     dxEffectiveRDFCSW_dVolume_dDensity(iComponent,jComponent) ) ) * dcEffectiveRDFCSW_dVolume(iComponent,jComponent) ) / &
      &     ( cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) ) )
    END IF
    ! Second derivative of the mean-attractive energy fluctuation with respect to the volume [K² . mol² / m⁶]
    IF( ZhangCorrectionLogical ) THEN
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * ZhangFactor ) + ( 2.D0 * dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * &
        &     ZhangFactor * dHSIsothermalCompressibility_dVolume ) + ( 2.D0 * rNumberDensity * dZhangFactor_dVolume * &
        &     dHSIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( 2.D0 * &
        &     dZhangFactor_dVolume * HSIsothermalCompressibility * dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) ) &
        &     + ( 2.D0 * HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * dZhangFactor_dVolume ) + ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     ZhangFactor * HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * &
        &     d2ZhangFactor_d2Volume * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * HSIsothermalCompressibility &
        &     ) + ( rNumberDensity * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
        &     d2HSIsothermalCompressibility_d2Volume * ZhangFactor ) + ( HSIsothermalCompressibility * &
        &     d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) * ZhangFactor ) + ( rNumberDensity * &
        &     HSIsothermalCompressibility * ZhangFactor * d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * ZhangFactor ) + ( 2.D0 * dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
        &     ZhangFactor * dHSIsothermalCompressibility_dVolume ) + ( 2.D0 * rNumberDensity * dZhangFactor_dVolume * &
        &     dHSIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( 2.D0 * &
        &     dZhangFactor_dVolume * HSIsothermalCompressibility * &
        &     dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) ) + ( 2.D0 * HSIsothermalCompressibility * &
        &     dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity * dZhangFactor_dVolume ) + &
        &     ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * ZhangFactor * HSIsothermalCompressibility * &
        &     dNumberDensity_dVolume ) + ( rNumberDensity * d2ZhangFactor_d2Volume * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * HSIsothermalCompressibility ) + ( rNumberDensity * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * d2HSIsothermalCompressibility_d2Volume * ZhangFactor &
        &     ) + ( HSIsothermalCompressibility * d2cMeanAttractiveEnergySL_d2Volume(iComponent,jComponent) * ZhangFactor ) &
        &     + ( rNumberDensity * HSIsothermalCompressibility * ZhangFactor * &
        &     d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * ZhangFactor ) + ( 2.D0 * dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * &
        &     ZhangFactor * dHSIsothermalCompressibility_dVolume ) + ( 2.D0 * rNumberDensity * dZhangFactor_dVolume * &
        &     dHSIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( 2.D0 * &
        &     dZhangFactor_dVolume * HSIsothermalCompressibility * &
        &     dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) ) + ( 2.D0 * HSIsothermalCompressibility * &
        &     dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * rNumberDensity * dZhangFactor_dVolume ) + &
        &     ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * ZhangFactor * HSIsothermalCompressibility * &
        &     dNumberDensity_dVolume ) + ( rNumberDensity * d2ZhangFactor_d2Volume * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * HSIsothermalCompressibility ) + ( rNumberDensity * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * d2HSIsothermalCompressibility_d2Volume * ZhangFactor &
        &     ) + ( HSIsothermalCompressibility * d2cMeanAttractiveEnergyYK_d2Volume(iComponent,jComponent) * ZhangFactor ) &
        &     + ( rNumberDensity * HSIsothermalCompressibility * ZhangFactor * &
        &     d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     rNumberDensity * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     dHCBIsothermalCompressibility_dVolume * ZhangFactor ) + ( 2.D0 * &
        &     dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * dHCBIsothermalCompressibility_dVolume * ZhangFactor ) + &
        &     ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * d2HCBIsothermalCompressibility_d2Volume * rNumberDensity * &
        &     ZhangFactor ) + ( d2cMeanAttractiveEnergyCSW_d2Volume(iComponent,jComponent) * HCBIsothermalCompressibility * &
        &     ZhangFactor ) + ( HCBIsothermalCompressibility * d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity * ZhangFactor ) + ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     HCBIsothermalCompressibility * dNumberDensity_dVolume * ZhangFactor ) + ( 2.D0 * rNumberDensity * &
        &     dHCBIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * dZhangFactor_dVolume &
        &     ) + ( 2.D0 * dZhangFactor_dVolume * dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * &
        &     HCBIsothermalCompressibility ) + ( 2.D0 * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     HCBIsothermalCompressibility * rNumberDensity * dZhangFactor_dVolume ) + ( rNumberDensity * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * d2ZhangFactor_d2Volume * HCBIsothermalCompressibility ) )
      END IF
    ELSE
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity ) + ( 2.D0 * dcMeanAttractiveEnergy_dVolume(iComponent,jComponent) * &
        &     dHSIsothermalCompressibility_dVolume ) + ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * d2HSIsothermalCompressibility_d2Volume * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( d2cMeanAttractiveEnergy_d2Volume(iComponent,jComponent) &
        &     * HSIsothermalCompressibility ) + ( rNumberDensity * d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) * &
        &     HSIsothermalCompressibility ) )
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity ) + ( 2.D0 * dcMeanAttractiveEnergySutherland_dVolume(iComponent,jComponent) * &
        &     dHSIsothermalCompressibility_dVolume ) + ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * d2HSIsothermalCompressibility_d2Volume * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( &
        &     d2cMeanAttractiveEnergySL_d2Volume(iComponent,jComponent) * HSIsothermalCompressibility ) + ( rNumberDensity * &
        &     d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) * HSIsothermalCompressibility ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     rNumberDensity ) + ( 2.D0 * dcMeanAttractiveEnergyYukawa_dVolume(iComponent,jComponent) * &
        &     dHSIsothermalCompressibility_dVolume ) + ( dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * d2HSIsothermalCompressibility_d2Volume * &
        &     dMeanAttractiveEnergy_dDensity(iComponent,jComponent) ) + ( &
        &     d2cMeanAttractiveEnergyYK_d2Volume(iComponent,jComponent) * HSIsothermalCompressibility ) + ( rNumberDensity * &
        &     d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) * HSIsothermalCompressibility ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        d2cMeanEnergyFluctuations_d2Volume(iComponent,jComponent) = 0.5D0 * ijaWellDepth(iComponent,jComponent) * ( ( 2.D0 * &
        &     rNumberDensity * dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * &
        &     dHCBIsothermalCompressibility_dVolume ) + ( 2.D0 * dcMeanAttractiveEnergyCSW_dVolume(iComponent,jComponent) * &
        &     dHCBIsothermalCompressibility_dVolume ) + ( dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
        &     d2HCBIsothermalCompressibility_d2Volume * rNumberDensity ) + ( HCBIsothermalCompressibility * &
        &     d2cMeanAttractiveEnergyCSW_d2Volume(iComponent,jComponent) ) + ( HCBIsothermalCompressibility * rNumberDensity * &
        &     d2FirstTPTCoeff_d2Volume_dDensity(iComponent,jComponent) ) + ( HCBIsothermalCompressibility * &
        &     dxFirstTPTCoefficient_dVolume_dDensity(iComponent,jComponent) * dNumberDensity_dVolume ) )
      END IF
    END IF
  END DO
END DO

! Mixing rule for the second-order perturbation coefficient [K² . mol² / m⁶]
CALL Mixing_Rules( mFraction, d2cMeanEnergyFluctuations_d2Volume, d2pSecondOrderCoefficient_d2Volume )

! Second derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MeanFluctuationFreeEnergy_d2Volume = d2pSecondOrderCoefficient_d2Volume * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  d2HigherOrderTerms_d2Volume = 0.D0
  DO nOrder = 3, nHigherOrder
    d2HigherOrderTerms_d2Volume = d2HigherOrderTerms_d2Volume + ( d2MeanAttractiveFreeEnergy_d2Volume * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) ) + ( dMeanAttractiveFreeEnergy_dVolume * ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * &
    &     ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy &
    &     * dMeanAttractiveFreeEnergy_dVolume) ) / ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) &
    &     + ( ( 4.D0 * DBLE( nOrder - 1 ) * DBLE( nOrder - 2 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 3 ) ) * &
    &     ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - &
    &     (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) * ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) ) / &
    &     ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) + &
    &     ( ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / &
    &     mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     d2MeanFluctuationFreeEnergy_d2Volume) - (mMeanAttFluctuationHelmholtzFreeEnergy * d2MeanAttractiveFreeEnergy_d2Volume) ) &
    &     / ( mMeanAttractiveHelmholtzFreeEnergy ) ) - ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) * &
    &     dMeanAttractiveFreeEnergy_dVolume ) / ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) )
  END DO
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
  d2HigherOrderFEnergy_d2Volume = d2HigherOrderTerms_d2Volume ! Proven units
ELSE
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
  d2HigherOrderFEnergy_d2Volume = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Second derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MonomerFreeEnergy_d2Volume = d2ExcludedVolumeFreeEnergy_d2Volume + d2MeanAttractiveFreeEnergy_d2Volume + &
&     d2MeanFluctuationFreeEnergy_d2Volume + d2HigherOrderFEnergy_d2Volume ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Second derivative of the total Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2TotalFreeEnergy_d2Volume = d2IdealFreeEnergy_d2Volume + d2MonomerFreeEnergy_d2Volume ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!  CROSS DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE AND THE VOLUME   !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Cross derivative of the ideal contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxIdealFreeEnergy_dVolume_dTemperature = dIdealFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Cross derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxEVFreeEnergy_dVolume_dTemperature = dExcludedVolumeFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Cross derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxSecondOrderTPT_dVolume_dTemperature = - dMeanFluctuationFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
IF( HigherOrderTPTLogical ) THEN
  ! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
  dMeanFluctuationFreeEnergy_dTemperature = - pSecondOrderCoefficient * cUniversalGas / ( Temperature * Temperature ) ! Proven units
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  dxHigherOrderTerms_dVolume_dTemperature = 0.D0
  DO nOrder = 3, nHigherOrder
    dxHigherOrderTerms_dVolume_dTemperature = dxHigherOrderTerms_dVolume_dTemperature + ( 2.D0 * dMeanAttractiveFreeEnergy_dVolume &
    &     * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / &
    &     mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * (dMeanFluctuationFreeEnergy_dTemperature / &
    &     mMeanAttractiveHelmholtzFreeEnergy) ) + ( 4.D0 * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 3 ) ) / &
    &     ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) * ( ( DBLE( nOrder - 2 ) * &
    &     dMeanFluctuationFreeEnergy_dTemperature * ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - &
    &     (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) ) + &
    &     ( mMeanAttFluctuationHelmholtzFreeEnergy * ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dxSecondOrderTPT_dVolume_dTemperature) - (dMeanFluctuationFreeEnergy_dTemperature * &
    &     dMeanAttractiveFreeEnergy_dVolume) ) ) )
  END DO
  ! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
  dxHigherOrderFEnergy_dVolume_dTemperature = dxHigherOrderTerms_dVolume_dTemperature ! Proven units
ELSE
  ! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
  dxHigherOrderFEnergy_dVolume_dTemperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Cross derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxMonomerFEnergy_dVolume_dTemperature = dxEVFreeEnergy_dVolume_dTemperature + dxSecondOrderTPT_dVolume_dTemperature + &
&     dxHigherOrderFEnergy_dVolume_dTemperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Cross derivative of the total Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxTotalFreeEnergy_dVolume_dTemperature = dxIdealFreeEnergy_dVolume_dTemperature + dxMonomerFEnergy_dVolume_dTemperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                     THERMODYNAMIC PROPERTIES                                     !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Isothermal compressibility [1 / Pa]
IsothermalCompressibility = 1.D0 / ( mVolume * d2TotalFreeEnergy_d2Volume )

! Thermal expansion coefficient [1 / K]
ThermalExpansionCoefficient = - IsothermalCompressibility * dxTotalFreeEnergy_dVolume_dTemperature

! Pressure [Pa]
Pressure = - dTotalFreeEnergy_dVolume

! Compressibility factor [unitless]
CompressibilityFactor = Pressure * mVolume / ( cUniversalGas * Temperature )

RETURN

END SUBROUTINE Calculate_Pressure

! ************************************************************************************************ !
! Calculates pressure and the derivatives of the Helmholtz free energy with respect to the volume  !
! for single components                                                                            !
! ************************************************************************************************ !
SUBROUTINE Calculate_Pressure_Single_Component( cComponent, mVolume, Temperature, Pressure, IsothermalCompressibility, &
&                                               ThermalExpansionCoefficient, CompressibilityFactor )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index
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
REAL( Kind= Real64 )                    :: Pressure                                  ! Pressure
REAL( Kind= Real64 )                    :: IsothermalCompressibility                 ! Isothermal compressibility
REAL( Kind= Real64 )                    :: ThermalExpansionCoefficient               ! Thermal expansion coefficient
REAL( Kind= Real64 )                    :: CompressibilityFactor                     ! Compressibility factor
REAL( Kind= Real64 )                    :: mVolume                                   ! Molar volume
REAL( Kind= Real64 )                    :: cMolecularVolume                          ! Molecular volume (component) [Å³]
REAL( Kind= Real64 )                    :: cSurfaceArea                              ! Surface area (component) [Å²]
REAL( Kind= Real64 )                    :: cCurvatureRadius                          ! Mean radius of curvature (component) [Å]
REAL( Kind= Real64 )                    :: cRadius                                   ! Radius (component) [Å]
REAL( Kind= Real64 )                    :: cNonSphericity                            ! Nonsphericity parameter (component)
REAL( Kind= Real64 )                    :: cSecondVirialCoefficient                  ! Second virial coefficient of the component
REAL( Kind= Real64 )                    :: fMolecularVolume                          ! Molecular volume (potential) [Å³]
REAL( Kind= Real64 )                    :: fSurfaceArea                              ! Surface area (potential) [Å²]
REAL( Kind= Real64 )                    :: fCurvatureRadius                          ! Mean radius of curvature (potential) [Å]
REAL( Kind= Real64 )                    :: fRadius                                   ! Radius (potential) [Å]
REAL( Kind= Real64 )                    :: fNonSphericity                            ! Nonsphericity parameter (potential)
REAL( Kind= Real64 )                    :: fSecondVirialCoefficient                  ! Second virial coefficient of the field
REAL( Kind= Real64 )                    :: MomentInertia                             ! Moment of inertia (non-spherical rigid bodies)
REAL( Kind= Real64 )                    :: Temperature                               ! Temperature
REAL( Kind= Real64 )                    :: rNumberDensity                            ! Reduced number density
REAL( Kind= Real64 )                    :: dNumberDensity_dVolume                    ! First derivative of the number density with respect to the volume
REAL( Kind= Real64 )                    :: dInverseNumberDensity_dVolume             ! First derivative of the inverse number density with respect to the volume
REAL( Kind= Real64 )                    :: d2NumberDensity_d2Volume                  ! Second derivative of the number density with respect to the volume
REAL( Kind= Real64 )                    :: mIdealHelmholtzFreeEnergy                 ! Ideal contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mExcludedVolumeHelmholtzFreeEnergy        ! Excluded-volume contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMeanAttractiveHelmholtzFreeEnergy        ! Mean-attractive energy contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMeanAttFluctuationHelmholtzFreeEnergy    ! Mean-attractive energy fluctuation contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMonomerMonomerHelmholtzFreeEnergy        ! Monomer-monomer contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mTotalHelmholtzFreeEnergy                 ! Total Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: HSBoublikMansoori                         ! Boublik-Mansoori hard-sphere mixture term
REAL( Kind= Real64 )                    :: AngleAverage                              ! Angle average of the excluded volume of a pair of non-spherical rigid bodies
REAL( Kind= Real64 )                    :: ZhangFactor                               ! Zhang's correction factor
REAL( Kind= Real64 )                    :: AuxDiameterRelationship                   ! Auxiliary factor (diameter relationship)
REAL( Kind= Real64 )                    :: HSIsothermalCompressibility               ! Hard-sphere isothermal compressibility (Percus-Yevick expression)
REAL( Kind= Real64 )                    :: HCBIsothermalCompressibility              ! Hard convex-body isothermal compressibility (Boublik expression)
REAL( Kind= Real64 )                    :: dMeanFluctuationFreeEnergy_dTemperature   ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dIdealFreeEnergy_dVolume                  ! First derivative of the ideal contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dHSBoublikMansoori_dVolume                ! First derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume
REAL( Kind= Real64 )                    :: dExcludedVolumeFreeEnergy_dVolume         ! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dHSIsothermalCompressibility_dVolume      ! First derivative of the hard-sphere isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                    :: dHCBIsothermalCompressibility_dVolume     ! First derivative of the hard convex-body isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                    :: dMeanAttractiveFreeEnergy_dVolume         ! First derivative of the mean-attractive energy contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dMeanFluctuationFreeEnergy_dVolume        ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dMonomerFreeEnergy_dVolume                ! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dTotalFreeEnergy_dVolume                  ! First derivative of the total Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2IdealFreeEnergy_d2Volume                ! Second derivative of the ideal contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2ExcludedVolumeFreeEnergy_d2Volume       ! Second derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2HSBoublikMansoori_d2Volume              ! Second derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume
REAL( Kind= Real64 )                    :: d2MeanAttractiveFreeEnergy_d2Volume       ! Second derivative of the mean-attractive energy contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2HSIsothermalCompressibility_d2Volume    ! Second derivative of the hard-sphere isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                    :: d2HCBIsothermalCompressibility_d2Volume   ! Second derivative of the hard convex-body isothermal compressibility with respect to the volume
REAL( Kind= Real64 )                    :: d2MeanFluctuationFreeEnergy_d2Volume      ! Second derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2MonomerFreeEnergy_d2Volume              ! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2TotalFreeEnergy_d2Volume                ! Second derivative of the total Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: dxIdealFreeEnergy_dVolume_dTemperature    ! Cross derivative of the ideal contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxEVFreeEnergy_dVolume_dTemperature       ! Cross derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxSecondOrderTPT_dVolume_dTemperature     ! Cross derivative of the second-order perturbation coefficient with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxMonomerFEnergy_dVolume_dTemperature     ! Cross derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxTotalFreeEnergy_dVolume_dTemperature    ! Cross derivative of the total Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxHigherOrderTerms_dVolume_dTemperature   ! Cross derivative of the higher-order terms with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: dxHigherOrderFEnergy_dVolume_dTemperature ! Cross derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume and the temperature
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Translational        ! Thermal de Broglie wavelength
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Translational_Cb     ! Thermal de Broglie wavelength (cubic)
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Rotational           ! Thermal de Broglie wavelength (rotational contribution)
REAL( Kind= Real64 )                    :: HigherOrderTerms                          ! Higher-order perturbation terms
REAL( Kind= Real64 )                    :: Factorial                                 ! Factorial function
REAL( Kind= Real64 )                    :: mHigherOrderHelmholtzFreeEnergy           ! Higher-order contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunct_dDensity      ! First derivative of the radial distribution function with respect to the density
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctYK_dDensity    ! First derivative of the radial distribution function with respect to the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctSL_dDensity    ! First derivative of the radial distribution function with respect to the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctCSW_dDensity   ! First derivative of the radial distribution function with respect to the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: dMeanAttractiveEnergy_dDensity            ! First derivative of the mean-attractive energy with respect to the density
REAL( Kind= Real64 )                    :: cAverageDiameterRelationship              ! Diameter relationship between components (the product of diameters divided by the sum of diameters)
REAL( Kind= Real64 )                    :: cRadialDistributionFunction               ! Contact radial distribution function
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunction      ! Contact radial distribution function for an effective packing fraction
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionYK    ! Contact radial distribution function for an effective packing fraction (Yukawa potential)
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionSL    ! Contact radial distribution function for an effective packing fraction (Sutherland potential)
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionCSW   ! Contact radial distribution function for an effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergy                     ! Mean-attractive energy between components
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergySutherland           ! Mean-attractive energy between components (Sutherland potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyYukawa               ! Mean-attractive energy between components (Yukawa potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyCSW                  ! Mean-attractive energy between components (Convex square-well potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyFluctuations         ! Mean-attractive energy fluctuations between components
REAL( Kind= Real64 )                    :: EffPackingFraction                        ! Effective packing fraction
REAL( Kind= Real64 )                    :: EffPackingFractionYukawa                  ! Effective packing fraction (Yukawa potential)
REAL( Kind= Real64 )                    :: EffPackingFractionSutherland              ! Effective packing fraction (Sutherland potential)
REAL( Kind= Real64 )                    :: EffPackingFractionCSW                     ! Effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 )                    :: dEffPackingFraction_dDensity              ! First derivative of the effective packing fraction with respect to the density
REAL( Kind= Real64 )                    :: dEffPackingFractionYukawa_dDensity        ! First derivative of the effective packing fraction with respect to the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionSutherland_dDensity    ! First derivative of the effective packing fraction with respect to the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dEffPackingFraction_dVolume               ! First derivative of the effective packing fraction with respect to the volume
REAL( Kind= Real64 )                    :: dEffPackingFractionYK_dVolume             ! First derivative of the effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionCSW_dVolume            ! First derivative of the effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionSL_dVolume             ! First derivative of the effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionCSW_dDensity           ! First derivative of the effective packing fraction with respect to the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: dxEffPFraction_dVolume_dDensity           ! Cross derivative of the effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 )                    :: dxEffPFractionYK_dVolume_dDensity         ! Cross derivative of the effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dxEffPFractionSL_dVolume_dDensity         ! Cross derivative of the effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dxEffPFractionCSW_dVolume_dDensity        ! Cross derivative of the effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: dxEffectiveRDF_dVolume_dDensity           ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 )                    :: dxEffectiveRDFYukawa_dVolume_dDensity     ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dxEffectiveRDFSutherland_dVolume_dDensity ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dxEffectiveRDFCSW_dVolume_dDensity        ! Cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: dcEffectiveRDF_dVolume                    ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume
REAL( Kind= Real64 )                    :: dcEffectiveRDFYukawa_dVolume              ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: dcEffectiveRDFSutherland_dVolume          ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: dcEffectiveRDFCSW_dVolume                 ! First derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: dcMeanAttractiveEnergy_dVolume            ! First derivative of the mean-attractive energy with respect to the volume
REAL( Kind= Real64 )                    :: dcMeanAttractiveEnergySutherland_dVolume  ! First derivative of the mean-attractive energy with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: dcMeanAttractiveEnergyYukawa_dVolume      ! First derivative of the mean-attractive energy with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: dcMeanAttractiveEnergyCSW_dVolume         ! First derivative of the mean-attractive energy with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: dcMeanEnergyFluctuations_dVolume          ! First derivative of the mean-attractive energy fluctuations with respect to the volume
REAL( Kind= Real64 )                    :: dxFirstTPTCoefficient_dVolume_dDensity    ! Cross derivative of the first-order perturbation coefficient with respect to the volume and the density
REAL( Kind= Real64 )                    :: dcRDFunction_dVolume                      ! First derivative of the contact radial distribution function with respect to the volume
REAL( Kind= Real64 )                    :: dZhangFactor_dVolume                      ! First derivative of the Zhang's correction factor with respect to the volume
REAL( Kind= Real64 )                    :: gAux1, gAux2                              ! Auxiliars (nonsphericity expressions)
REAL( Kind= Real64 )                    :: dHigherOrderTerms_dVolume                 ! First derivative of the higher-order terms with respect to the volume
REAL( Kind= Real64 )                    :: d2HigherOrderTerms_d2Volume               ! Second derivative of the higher-order terms with respect to the volume
REAL( Kind= Real64 )                    :: dHigherOrderFEnergy_dVolume               ! First derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2HigherOrderFEnergy_d2Volume             ! Second derivative of the higher-order contribution to the Helmholtz free energy with respect to the volume
REAL( Kind= Real64 )                    :: d2cMeanAttractiveEnergy_d2Volume          ! Second derivative of the mean-attractive energy with respect to the volume
REAL( Kind= Real64 )                    :: d2cMeanAttractiveEnergySL_d2Volume        ! Second derivative of the mean-attractive energy with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: d2cMeanAttractiveEnergyYK_d2Volume        ! Second derivative of the mean-attractive energy with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: d2cMeanAttractiveEnergyCSW_d2Volume       ! Second derivative of the mean-attractive energy with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: d2cEffectiveRDF_d2Volume                  ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume
REAL( Kind= Real64 )                    :: d2cEffectiveRDFYukawa_d2Volume            ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: d2cEffectiveRDFSutherland_d2Volume        ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: d2cEffectiveRDFCSW_d2Volume               ! Second derivative of the contact radial distribution function for an effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: d2EffPackingFraction_d2Volume             ! Second derivative of the effective packing fraction with respect to the volume
REAL( Kind= Real64 )                    :: d2EffPackingFractionYK_d2Volume           ! Second derivative of the effective packing fraction with respect to the volume (Yukawa potential)
REAL( Kind= Real64 )                    :: d2EffPackingFractionSL_d2Volume           ! Second derivative of the effective packing fraction with respect to the volume (Sutherland potential)
REAL( Kind= Real64 )                    :: d2EffPackingFractionCSW_d2Volume          ! Second derivative of the effective packing fraction with respect to the volume (Convex square-well potential)
REAL( Kind= Real64 )                    :: d2EffPFraction_d2Volume_dDensity          ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 )                    :: d2EffPFractionYK_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 )                    :: d2EffPFractionSL_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 )                    :: d2EffPFractionCSW_d2Volume_dDensity       ! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: d2cRDFunction_d2Volume                    ! Second derivative of the contact radial distribution function with respect to the volume
REAL( Kind= Real64 )                    :: d2cMeanEnergyFluctuations_d2Volume        ! Second derivative of the mean-attractive energy fluctuations with respect to the volume
REAL( Kind= Real64 )                    :: d2FirstTPTCoeff_d2Volume_dDensity         ! First derivative with respect to the volume of the cross derivative of the first-order perturbation coefficient with respect to the volume and the density
REAL( Kind= Real64 )                    :: d2EffectiveRDF_d2Volume_dDensity          ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 )                    :: d2EffectiveRDFYK_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 )                    :: d2EffectiveRDFSL_d2Volume_dDensity        ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 )                    :: d2EffectiveRDFCSW_d2Volume_dDensity       ! First derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: d2ZhangFactor_d2Volume                    ! Second derivative of the Zhang's correction factor with respect to the volume
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxEffRDF_d2Volume_dDensity             ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxEffRDFYukawa_d2Volume_dDensity       ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxEffRDFSutherland_d2Volume_dDensity   ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxEffRDFCSW_d2Volume_dDensity          ! Auxiliary factor of a first derivative with respect to the volume of the cross derivative of the contact radial distribution function for an effective packing fraction with respect to the volume and the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficients            ! Coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsYukawa      ! Coefficients of the effective packing fraction (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsSutherland  ! Coefficients of the effective packing fraction (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsCSW         ! Coefficients of the effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: dAuxiliaryHSFactor_dVolume                ! Auxiliary factor of a first derivative with respect to the volume (hard-sphere contribution)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxiliaryHSFactor_d2Volume              ! Auxiliary factor of a second derivative with respect to the volume (hard-sphere contribution)
REAL( Kind= Real64 ), DIMENSION( 6 )    :: d2AuxHSIC_d2Volume                        ! Auxiliary factor of a second derivative with respect to the volume (hard-sphere isothermal compressibility)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: d2AuxHCBIC_d2Volume                       ! Auxiliary factor of a second derivative with respect to the volume (hard convex-body isothermal compressibility)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: drDensity_dVolume                         ! First derivative of the reduced densities with respect to the volume
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: d2rDensity_d2Volume                       ! Second derivative of the reduced densities with respect to the volume
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: rDensity                                  ! Reduced densities
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: rDensityConstants                         ! Reduced density constants
REAL( Kind= Real64 ), DIMENSION( 2, 9 ) :: CSWAlphaCoefficients                      ! Coefficients of the nonsphericity (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ) :: CSW1CoefficientsMatrix                    ! Matrix of coefficients of the effective packing fraction (first-order coefficient)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ) :: CSW2CoefficientsMatrix                    ! Matrix of coefficients of the effective packing fraction (second-order coefficient)

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
mExcludedVolumeHelmholtzFreeEnergy = HSBoublikMansoori * AngleAverage * cUniversalGas * Temperature ! Proven units
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
  EffPackingFractionSutherland = EffPackingFractionCoefficientsSutherland(1) * rDensity(3) + &
  &     EffPackingFractionCoefficientsSutherland(2) * rDensity(3) * rDensity(3)
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
  &     EffPackingFraction * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * &
  &     ( 1.D0 - EffPackingFraction ) ) )
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
!            FIRST DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE VOLUME             !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! First derivative of the number density with respect to the volume [mol / (Å³ . m³)]
dNumberDensity_dVolume = - rNumberDensity / mVolume

! First derivative of the inverse number density with respect to the volume [mol . (Å³ / m³)]
dInverseNumberDensity_dVolume = 1.D0 / rNumberDensity / mVolume

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dIdealFreeEnergy_dVolume = - (1.D0 / mVolume) * cUniversalGas * Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! First derivative of the reduced densities with respect to the volume
drDensity_dVolume(0) = rDensityConstants(0) * dNumberDensity_dVolume ! [mol / (Å³ . m³)]
drDensity_dVolume(1) = rDensityConstants(1) * dNumberDensity_dVolume ! [mol / (Å² . m³)]
drDensity_dVolume(2) = rDensityConstants(2) * dNumberDensity_dVolume ! [mol / (Å . m³)]
drDensity_dVolume(3) = rDensityConstants(3) * dNumberDensity_dVolume ! [mol / m³]

! Auxiliary factors of the derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume ! [mol / (Å³ . m³)]
dAuxiliaryHSFactor_dVolume(1) = - ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * &
&     ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * drDensity_dVolume(3) + DLOG( 1.D0 - rDensity(3) ) * ( ( 3.D0 * rDensity(2) * rDensity(2) &
&     * drDensity_dVolume(2) / ( rDensity(3) * rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * &
&     drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) * rDensity(3) ) ) - drDensity_dVolume(0) )
dAuxiliaryHSFactor_dVolume(2) = 3.D0 * ( drDensity_dVolume(1) * rDensity(2) + drDensity_dVolume(2) * rDensity(1) ) / ( 1.D0 - &
&     rDensity(3) ) + 3.D0 * rDensity(1) * rDensity(2) * drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) )
dAuxiliaryHSFactor_dVolume(3) = ( rDensity(2) * rDensity(2) * ( 3.D0 * drDensity_dVolume(2) * rDensity(3) * rDensity(3) - &
&     ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * drDensity_dVolume(2) ) * rDensity(3) + drDensity_dVolume(3) * &
&     rDensity(2) ) ) / ( ( rDensity(3) - 1.D0 ) * ( rDensity(3) - 1.D0 ) * ( rDensity(3) - 1.D0 ) * rDensity(3) * rDensity(3) )

! First derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume [mol / m³]
dHSBoublikMansoori_dVolume = HSBoublikMansoori * rNumberDensity * dInverseNumberDensity_dVolume + ( 6.D0 / (cPi * &
&     rNumberDensity) ) * SUM( dAuxiliaryHSFactor_dVolume )

! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dExcludedVolumeFreeEnergy_dVolume = dHSBoublikMansoori_dVolume * AngleAverage * cUniversalGas * Temperature ! Proven units
IF( ReferenceBoublikLogical ) THEN
  dExcludedVolumeFreeEnergy_dVolume = - ( ( ( 1.D0 + 3.D0 * cNonSphericity ) * rDensity(3) ) + ( ( 3.D0 * cNonSphericity * &
  &     cNonSphericity - 3.D0 * cNonSphericity - 2.D0 ) * rDensity(3) * rDensity(3) ) + ( ( 1.D0 - cNonSphericity * &
  &     cNonSphericity ) * rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) )
  dExcludedVolumeFreeEnergy_dVolume = dExcludedVolumeFreeEnergy_dVolume * cUniversalGas * Temperature
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! First derivative of the effective packing fraction with respect to the volume [mol / m³]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dEffPackingFraction_dVolume = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dEffPackingFraction_dVolume = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3)
  dEffPackingFractionSL_dVolume = EffPackingFractionCoefficientsSutherland(1) * rDensity(3) + 2.D0 * &
  &     EffPackingFractionCoefficientsSutherland(2) * rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dEffPackingFraction_dVolume = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3)
  dEffPackingFractionYK_dVolume = EffPackingFractionCoefficientsYukawa(1) * rDensity(3) + 2.D0 * &
  &     EffPackingFractionCoefficientsYukawa(2) * rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dEffPackingFraction_dVolume = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3)
  dEffPackingFractionCSW_dVolume = EffPackingFractionCoefficientsCSW(1) * rDensity(3) + 2.D0 * rDensity(3) * rDensity(3) * &
  &     EffPackingFractionCoefficientsCSW(2) + 3.D0 * rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsCSW(3)
END IF
dEffPackingFraction_dVolume = - dEffPackingFraction_dVolume / mVolume
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dEffPackingFractionSL_dVolume = - dEffPackingFractionSL_dVolume / mVolume
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dEffPackingFractionYK_dVolume = - dEffPackingFractionYK_dVolume / mVolume
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dEffPackingFractionCSW_dVolume = - dEffPackingFractionCSW_dVolume / mVolume
END IF

! Cross derivative of the effective packing fraction with respect to the volume and the density [mol . (Å³ / m³)]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dxEffPFraction_dVolume_dDensity = dEffPackingFraction_dDensity / mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * &
  &     ( ( EffPackingFractionCoefficients(1) * rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3) ) + ( 9.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3) * rDensity(3) ) )
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dxEffPFraction_dVolume_dDensity = dEffPackingFraction_dDensity / mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * &
  &     ( ( EffPackingFractionCoefficients(1) * rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3) ) )
  dxEffPFractionSL_dVolume_dDensity = dEffPackingFractionSutherland_dDensity / mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * &
  &     ( ( EffPackingFractionCoefficientsSutherland(1) * rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficientsSutherland(2) * &
  &     rDensity(3) * rDensity(3) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dxEffPFraction_dVolume_dDensity = dEffPackingFraction_dDensity / mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * &
  &     ( ( EffPackingFractionCoefficients(1) * rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3) ) )
  dxEffPFractionYK_dVolume_dDensity = dEffPackingFractionYukawa_dDensity / mVolume - ( 1.D0 / ( rNumberDensity * mVolume ) ) * ( ( &
  &     EffPackingFractionCoefficientsYukawa(1) * rDensity(3) ) + ( 4.D0 * EffPackingFractionCoefficientsYukawa(2) * rDensity(3) * &
  &     rDensity(3) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dxEffPFraction_dVolume_dDensity = - ( ( 2.D0 * EffPackingFractionCoefficients(2) * rDensity(3) * rDensity(3) ) + ( 6.D0 * &
  &     EffPackingFractionCoefficients(3) * rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * rNumberDensity )
  dxEffPFractionCSW_dVolume_dDensity = - ( ( 2.D0 * EffPackingFractionCoefficientsCSW(2) * rDensity(3) * rDensity(3) ) + ( 6.D0 * &
  &     EffPackingFractionCoefficientsCSW(3) * rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * rNumberDensity )
END IF

! First derivative of the contact radial distribution function with respect to the volume [mol / m³]
dcRDFunction_dVolume = ( drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) * ( 1.D0 + ( 3.D0 * &
&     AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + 2.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) + ( 2.D0 * &
&     AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
&     rDensity(3) / ( 1.D0 - rDensity(3) ) ) + ( 3.D0 * rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) ) ) ) )

! First derivative of the radial distribution function for an effective packing fraction with respect to the volume [mol / m³]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  dcEffectiveRDF_dVolume = ( dEffPackingFraction_dVolume / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) * &
  &     ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + 2.D0 * EffPackingFraction / ( 1.D0 - &
  &     EffPackingFraction ) ) ) + ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * &
  &     cAverageDiameterRelationship * ( ( 2.D0 * EffPackingFraction / ( 1.D0 - EffPackingFraction ) ) + ( 3.D0 * &
  &     EffPackingFraction * EffPackingFraction / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) ) ) )
ELSE
  dcEffectiveRDF_dVolume = ( dEffPackingFraction_dVolume / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) ) + ( 3.D0 * gAux1 * ( dEffPackingFraction_dVolume * ( 1.D0 + EffPackingFraction ) ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + ( ( 2.D0 * gAux2 * &
  &     EffPackingFraction * ( ( 2.D0 + EffPackingFraction ) * dEffPackingFraction_dVolume ) ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dcEffectiveRDFSutherland_dVolume = ( dEffPackingFractionSL_dVolume / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + &
  &     2.D0 * EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( 2.D0 * AuxDiameterRelationship * &
  &     AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
  &     EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) + ( 3.D0 * EffPackingFractionSutherland * &
  &     EffPackingFractionSutherland / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dcEffectiveRDFYukawa_dVolume = ( dEffPackingFractionYK_dVolume / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + 2.D0 &
  &     * EffPackingFractionYukawa / ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( 2.D0 * AuxDiameterRelationship * &
  &     AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
  &     EffPackingFractionYukawa / ( 1.D0 - EffPackingFractionYukawa ) ) + ( 3.D0 * EffPackingFractionYukawa * &
  &     EffPackingFractionYukawa / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dcEffectiveRDFCSW_dVolume = ( dEffPackingFractionCSW_dVolume * ( 1.D0 + 3.D0 * gAux1 ) / ( ( 1.D0 - EffPackingFractionCSW ) * &
  &     ( 1.D0 - EffPackingFractionCSW ) ) ) + ( 2.D0 * EffPackingFractionCSW * dEffPackingFractionCSW_dVolume * ( 3.D0 * gAux1 + &
  &     2.D0 * gAux2 ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) &
  &     ) ) + ( 6.D0 * EffPackingFractionCSW * EffPackingFractionCSW * gAux2 * dEffPackingFractionCSW_dVolume / ( ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) )
END IF

! First derivative of the mean-attractive energy with respect to the volume [K . mol / m³]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dcMeanAttractiveEnergy_dVolume = - 2.D0 * cPi * ( ijPotentialRangeCubic(cComponent,cComponent) - 1.D0 ) * &
  &     ( ijaDiameterSphereCubic(cComponent,cComponent) / 3.D0 ) * ( ( dNumberDensity_dVolume * &
  &     cEffectiveRadialDistributionFunction ) + dcEffectiveRDF_dVolume * rNumberDensity ) * ijaWellDepth(cComponent,cComponent)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dcMeanAttractiveEnergy_dVolume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) * &
  &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunction ) + ( rNumberDensity * dcEffectiveRDF_dVolume ) )
  dcMeanAttractiveEnergySutherland_dVolume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) * &
  &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunctionSL ) + ( rNumberDensity * &
  &     dcEffectiveRDFSutherland_dVolume ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dcMeanAttractiveEnergy_dVolume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) + ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) ) ) * ( ( dNumberDensity_dVolume * &
  &     cEffectiveRadialDistributionFunction ) + ( rNumberDensity * dcEffectiveRDF_dVolume ) )
  dcMeanAttractiveEnergyYukawa_dVolume = - cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) * &
  &     ( ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunctionYK ) + ( rNumberDensity * dcEffectiveRDFYukawa_dVolume ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dcMeanAttractiveEnergy_dVolume = - ijaWellDepth(cComponent,cComponent) * ( fMolecularVolume * ( 1.D0 + 3.D0 * fNonSphericity ) - &
  &     cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) ) * ( dNumberDensity_dVolume * cEffectiveRadialDistributionFunction + &
  &     rNumberDensity * dcEffectiveRDF_dVolume )
  dcMeanAttractiveEnergyCSW_dVolume = - ijaWellDepth(cComponent,cComponent) * ( fMolecularVolume * ( 1.D0 + 3.D0 * &
  &     fNonSphericity ) - cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) ) * ( dNumberDensity_dVolume * &
  &     cEffectiveRadialDistributionFunctionCSW + rNumberDensity * dcEffectiveRDFCSW_dVolume )
END IF

! First derivative of the first-order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMeanAttractiveFreeEnergy_dVolume = dcMeanAttractiveEnergy_dVolume * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the hard-sphere isothermal compressibility with respect to the volume [mol / m³]
dHSIsothermalCompressibility_dVolume = ( ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
&     * drDensity_dVolume(3) * rDensity(0) ) ) / ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * &
&     rDensity(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ) ) - &
&     ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) / &
&     ( ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ) * ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * &
&     rDensity(2) ) ) ) ) * ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * rDensity(0) * drDensity_dVolume(3) ) + ( 6.D0 * ( ( drDensity_dVolume(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( drDensity_dVolume(2) * rDensity(1) * ( 1.D0 - rDensity(3) ) ) - ( drDensity_dVolume(3) * rDensity(1) * &
&     rDensity(2) ) ) ) + ( 27.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) )

! First derivative of the hard convex body isothermal compressibility with respect to the volume [mol / m³]
IF( PYHCBCorrectionLogical ) THEN
  dHCBIsothermalCompressibility_dVolume = ( ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
  &     rDensity(3) ) / ( mVolume * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) ) ) ) + ( ( 2.D0 * rDensity(3) * ( 3.D0 * &
  &     cNonSphericity - 1.D0 ) + 2.D0 * rDensity(3) * rDensity(3) * ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * &
  &     cNonSphericity + 1.D0 ) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) ) / ( mVolume * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) ) * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * &
  &     cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + &
  &     1.D0 ) ) )
ELSE
  dHCBIsothermalCompressibility_dVolume = ( ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
  &     rDensity(3) ) / ( mVolume * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 4.D0 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) * cNonSphericity * cNonSphericity + rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity ) ) ) + ( ( 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + 2.D0 * rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 12.D0 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) * cNonSphericity * cNonSphericity + 4.D0 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     cNonSphericity * cNonSphericity ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) ) / ( mVolume * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 4.D0 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) * cNonSphericity * cNonSphericity + rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity ) * ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * &
  &     ( 9.D0 * cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 4.D0 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) * cNonSphericity * cNonSphericity + rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity ) )
END IF

! First derivative of the Zhang's factor with respect to the volume [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    dZhangFactor_dVolume = 4.D0 * ZhangCorrection * rDensity(3) * drDensity_dVolume(3)
  ELSE
    dZhangFactor_dVolume = - 0.25D0 * ZhangCorrection * ijSecondVirialCoefficient(cComponent,cComponent) * &
    &     ijSecondVirialCoefficient(cComponent,cComponent) * rNumberDensity * rNumberDensity / mVolume
  END IF
END IF

! Cross derivative of the effective radial distribution function with respect to the volume and the density [mol . (Å³ / m³)]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  dxEffectiveRDF_dVolume_dDensity = ( ( dxEffPFraction_dVolume_dDensity * ( 1.D0 - EffPackingFraction ) + ( 2.D0 * &
  &     dEffPackingFraction_dDensity * dEffPackingFraction_dVolume ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + ( 3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship &
  &     * ( ( ( ( dxEffPFraction_dVolume_dDensity * ( 1.D0 + EffPackingFraction ) ) + ( dEffPackingFraction_dDensity * &
  &     dEffPackingFraction_dVolume ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) ) + ( ( 3.D0 * dEffPackingFraction_dDensity * ( 1.D0 + EffPackingFraction ) * &
  &     dEffPackingFraction_dVolume ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) ) ) + ( 2.D0 * cAverageDiameterRelationship * &
  &     cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 * &
  &     dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * ( 1.D0 + EffPackingFraction ) ) + ( ( 2.D0 + &
  &     EffPackingFraction ) * EffPackingFraction * dxEffPFraction_dVolume_dDensity ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 &
  &     - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + ( ( 4.D0 * ( 2.D0 + &
  &     EffPackingFraction ) * EffPackingFraction * dEffPackingFraction_dVolume * dEffPackingFraction_dDensity ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * &
  &     ( 1.D0 - EffPackingFraction ) ) ) ) )
ELSE
  dxEffectiveRDF_dVolume_dDensity = ( dxEffPFraction_dVolume_dDensity / ( ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) ) ) + ( ( dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * ( 3.D0 * gAux1 + &
  &     2.D0 ) + 3.D0 * gAux1 * dxEffPFraction_dVolume_dDensity * ( 1.D0 + EffPackingFraction ) ) / &
  &     ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + &
  &     ( ( ( dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * ( ( 9.D0 * gAux1 + 4.D0 * gAux2 ) * ( 1.D0 + &
  &     EffPackingFraction ) ) ) + 2.D0 * gAux2 * EffPackingFraction * dxEffPFraction_dVolume_dDensity * ( 2.D0 + &
  &     EffPackingFraction ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) ) ) + ( ( 8.D0 * gAux2 * EffPackingFraction * ( 2.D0 + EffPackingFraction ) * &
  &     dEffPackingFraction_dVolume * dEffPackingFraction_dDensity ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dxEffectiveRDFSutherland_dVolume_dDensity = ( ( dxEffPFractionSL_dVolume_dDensity * ( 1.D0 - EffPackingFractionSutherland ) + ( &
  &     2.D0 * dEffPackingFractionSutherland_dDensity * dEffPackingFractionSL_dVolume ) ) / ( ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( &
  &     3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship * ( ( ( ( dxEffPFractionSL_dVolume_dDensity * ( 1.D0 + &
  &     EffPackingFractionSutherland ) ) + ( dEffPackingFractionSutherland_dDensity * dEffPackingFractionSL_dVolume ) ) / ( ( 1.D0 &
  &     - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) + &
  &     ( ( 3.D0 * dEffPackingFractionSutherland_dDensity * ( 1.D0 + EffPackingFractionSutherland ) * &
  &     dEffPackingFractionSL_dVolume ) / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) ) ) + ( 2.D0 * &
  &     cAverageDiameterRelationship * cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( &
  &     2.D0 * dEffPackingFractionSL_dVolume * dEffPackingFractionSutherland_dDensity * ( 1.D0 + EffPackingFractionSutherland ) ) &
  &     + ( ( 2.D0 + EffPackingFractionSutherland ) * EffPackingFractionSutherland * dxEffPFractionSL_dVolume_dDensity ) ) / ( ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) &
  &     * ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( ( 4.D0 * ( 2.D0 + EffPackingFractionSutherland ) * &
  &     EffPackingFractionSutherland * dEffPackingFractionSL_dVolume * dEffPackingFractionSutherland_dDensity ) / ( ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dxEffectiveRDFYukawa_dVolume_dDensity = ( ( dxEffPFractionYK_dVolume_dDensity * ( 1.D0 - EffPackingFractionYukawa ) + ( 2.D0 * &
  &     dEffPackingFractionYukawa_dDensity * dEffPackingFractionYK_dVolume ) ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( 3.D0 * cAverageDiameterRelationship * &
  &     AuxDiameterRelationship * ( ( ( ( dxEffPFractionYK_dVolume_dDensity * ( 1.D0 + EffPackingFractionYukawa ) ) + ( &
  &     dEffPackingFractionYukawa_dDensity * dEffPackingFractionYK_dVolume ) ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( ( 3.D0 * dEffPackingFractionYukawa_dDensity * ( &
  &     1.D0 + EffPackingFractionYukawa ) * dEffPackingFractionYK_dVolume ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) ) ) + ( 2.D0 * &
  &     cAverageDiameterRelationship * cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( &
  &     2.D0 * dEffPackingFractionYK_dVolume * dEffPackingFractionYukawa_dDensity * ( 1.D0 + EffPackingFractionYukawa ) ) + ( ( &
  &     2.D0 + EffPackingFractionYukawa ) * EffPackingFractionYukawa * dxEffPFractionYK_dVolume_dDensity ) ) / ( ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) + ( ( 4.D0 * ( 2.D0 + EffPackingFractionYukawa ) * EffPackingFractionYukawa * &
  &     dEffPackingFractionYK_dVolume * dEffPackingFractionYukawa_dDensity ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dxEffectiveRDFCSW_dVolume_dDensity = ( dxEffPFractionCSW_dVolume_dDensity / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) ) + ( ( dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dDensity * ( 3.D0 * gAux1 + &
  &     2.D0 ) + 3.D0 * gAux1 * dxEffPFractionCSW_dVolume_dDensity * ( 1.D0 + EffPackingFractionCSW ) ) / &
  &     ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) ) + &
  &     ( ( ( dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dDensity * ( ( 9.D0 * gAux1 + 4.D0 * gAux2 ) * ( 1.D0 + &
  &     EffPackingFractionCSW ) ) ) + 2.D0 * gAux2 * EffPackingFractionCSW * dxEffPFractionCSW_dVolume_dDensity * ( 2.D0 + &
  &     EffPackingFractionCSW ) ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) ) + ( ( 8.D0 * gAux2 * EffPackingFractionCSW * ( 2.D0 + &
  &     EffPackingFractionCSW ) * dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dDensity ) / ( ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) )
END IF

! Cross derivative of the first-order perturbation coefficient with respect to the volume and the density [K . mol . (Å³ / m³)]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dxFirstTPTCoefficient_dVolume_dDensity = ( dcMeanAttractiveEnergy_dVolume / rNumberDensity ) - ( cMeanAttractiveEnergy * &
  &     dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) + ( ( ( dcMeanAttractiveEnergy_dVolume * &
  &     dEffRadialDistributionFunct_dDensity ) - ( ( cMeanAttractiveEnergy / cEffectiveRadialDistributionFunction ) * &
  &     dcEffectiveRDF_dVolume * dEffRadialDistributionFunct_dDensity ) + ( cMeanAttractiveEnergy * &
  &     dxEffectiveRDF_dVolume_dDensity ) ) * ( 1.D0 / cEffectiveRadialDistributionFunction ) )
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dxFirstTPTCoefficient_dVolume_dDensity = ( dcMeanAttractiveEnergySutherland_dVolume / rNumberDensity ) - &
  &     ( cMeanAttractiveEnergySutherland * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) + ( ( ( &
  &     dcMeanAttractiveEnergySutherland_dVolume * dEffRadialDistributionFunctSL_dDensity ) - ( ( cMeanAttractiveEnergySutherland &
  &     / cEffectiveRadialDistributionFunctionSL ) * dcEffectiveRDFSutherland_dVolume * dEffRadialDistributionFunctSL_dDensity ) + &
  &     ( cMeanAttractiveEnergySutherland * dxEffectiveRDFSutherland_dVolume_dDensity ) ) * ( 1.D0 / &
  &     cEffectiveRadialDistributionFunctionSL ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dxFirstTPTCoefficient_dVolume_dDensity = ( dcMeanAttractiveEnergyYukawa_dVolume / rNumberDensity ) - &
  &     ( cMeanAttractiveEnergyYukawa * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) + ( ( ( &
  &     dcMeanAttractiveEnergyYukawa_dVolume * dEffRadialDistributionFunctYK_dDensity ) - ( ( cMeanAttractiveEnergyYukawa / &
  &     cEffectiveRadialDistributionFunctionYK ) * dcEffectiveRDFYukawa_dVolume * dEffRadialDistributionFunctYK_dDensity ) &
  &     + ( cMeanAttractiveEnergyYukawa * dxEffectiveRDFYukawa_dVolume_dDensity ) ) * ( 1.D0 / &
  &     cEffectiveRadialDistributionFunctionYK ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dxFirstTPTCoefficient_dVolume_dDensity = ( dcMeanAttractiveEnergyCSW_dVolume / rNumberDensity ) - ( cMeanAttractiveEnergyCSW * &
  &     dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) + ( ( ( dcMeanAttractiveEnergyCSW_dVolume * &
  &     dEffRadialDistributionFunctCSW_dDensity ) - ( ( cMeanAttractiveEnergyCSW / cEffectiveRadialDistributionFunctionCSW ) * &
  &     dcEffectiveRDFCSW_dVolume * dEffRadialDistributionFunctCSW_dDensity ) + ( cMeanAttractiveEnergyCSW * &
  &     dxEffectiveRDFCSW_dVolume_dDensity ) ) / cEffectiveRadialDistributionFunctionCSW )
END IF

! First derivative of the mean-attractive energy fluctuation with respect to the volume [K² . mol / m³]
IF( ZhangCorrectionLogical ) THEN
  IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity * ZhangFactor ) + ( dcMeanAttractiveEnergy_dVolume * &
    &     HSIsothermalCompressibility * ZhangFactor ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     ZhangFactor * rNumberDensity ) + ( HSIsothermalCompressibility * dMeanAttractiveEnergy_dDensity * rNumberDensity * &
    &     dZhangFactor_dVolume ) )
  ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity * ZhangFactor ) + ( dcMeanAttractiveEnergySutherland_dVolume * &
    &     HSIsothermalCompressibility * ZhangFactor ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     ZhangFactor * rNumberDensity ) + ( HSIsothermalCompressibility * dMeanAttractiveEnergy_dDensity * rNumberDensity * &
    &     dZhangFactor_dVolume ) )
  ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity * ZhangFactor ) + ( dcMeanAttractiveEnergyYukawa_dVolume * &
    &     HSIsothermalCompressibility * ZhangFactor ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     ZhangFactor * rNumberDensity ) + ( HSIsothermalCompressibility * dMeanAttractiveEnergy_dDensity * rNumberDensity * &
    &     dZhangFactor_dVolume ) )
  ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHCBIsothermalCompressibility_dVolume * rNumberDensity * ZhangFactor ) + ( dcMeanAttractiveEnergyCSW_dVolume * &
    &     HCBIsothermalCompressibility * ZhangFactor ) + ( HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     rNumberDensity * ZhangFactor ) + ( HCBIsothermalCompressibility * dMeanAttractiveEnergy_dDensity * rNumberDensity * &
    &     dZhangFactor_dVolume ) )
  END IF
ELSE
  IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity ) + ( dcMeanAttractiveEnergy_dVolume * HSIsothermalCompressibility &
    &     ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity ) )
  ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity ) + ( dcMeanAttractiveEnergySutherland_dVolume * &
    &     HSIsothermalCompressibility ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     rNumberDensity ) )
  ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHSIsothermalCompressibility_dVolume * rNumberDensity ) + ( dcMeanAttractiveEnergyYukawa_dVolume * &
    &     HSIsothermalCompressibility ) + ( HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     rNumberDensity ) )
  ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
    dcMeanEnergyFluctuations_dVolume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( dMeanAttractiveEnergy_dDensity * &
    &     dHCBIsothermalCompressibility_dVolume * rNumberDensity ) + ( dcMeanAttractiveEnergyCSW_dVolume * &
    &     HCBIsothermalCompressibility ) + ( HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     rNumberDensity ) )
  END IF
END IF

! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMeanFluctuationFreeEnergy_dVolume = dcMeanEnergyFluctuations_dVolume * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  dHigherOrderTerms_dVolume = 0.D0
  DO nOrder = 3, nHigherOrder
    dHigherOrderTerms_dVolume = dHigherOrderTerms_dVolume + ( dMeanAttractiveFreeEnergy_dVolume * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) ) + ( ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy &
    &     * dMeanAttractiveFreeEnergy_dVolume) ) / ( mMeanAttractiveHelmholtzFreeEnergy ) ) )
  END DO
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
  dHigherOrderFEnergy_dVolume = dHigherOrderTerms_dVolume ! Proven units
ELSE
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
  dHigherOrderFEnergy_dVolume = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! First derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dMonomerFreeEnergy_dVolume = dExcludedVolumeFreeEnergy_dVolume + dMeanAttractiveFreeEnergy_dVolume + &
&     dMeanFluctuationFreeEnergy_dVolume + dHigherOrderFEnergy_dVolume ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! First derivative of the total Helmholtz free energy with respect to the volume [Joule / m³ or Pa]
dTotalFreeEnergy_dVolume = dIdealFreeEnergy_dVolume + dMonomerFreeEnergy_dVolume ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!            SECOND DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE VOLUME            !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Second derivative of the number density with respect to the volume [mol² / (Å³ . m⁶)]
d2NumberDensity_d2Volume = 2.D0 * rNumberDensity / ( mVolume * mVolume )

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2IdealFreeEnergy_d2Volume = ( 1.D0 / ( mVolume * mVolume ) ) * cUniversalGas * Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Second derivative of the reduced densities with respect to the volume
d2rDensity_d2Volume(0) = rDensityConstants(0) * d2NumberDensity_d2Volume ! [mol² / (Å³ . m⁶)]
d2rDensity_d2Volume(1) = rDensityConstants(1) * d2NumberDensity_d2Volume ! [mol² / (Å² . m⁶)]
d2rDensity_d2Volume(2) = rDensityConstants(2) * d2NumberDensity_d2Volume ! [mol² / (Å . m⁶)]
d2rDensity_d2Volume(3) = rDensityConstants(3) * d2NumberDensity_d2Volume ! [mol² / m⁶]

! Auxiliary factors of the derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume ! [mol² / (Å³ . m⁶)]
d2AuxiliaryHSFactor_d2Volume(1) = - ( ( ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) &
&     * rDensity(3) ) ) - drDensity_dVolume(0) ) * ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * drDensity_dVolume(3) ) - ( ( ( rDensity(2) &
&     * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * drDensity_dVolume(3) * drDensity_dVolume(3) &
&     / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * &
&     rDensity(3) ) ) - rDensity(0) ) * ( 1.D0 / ( 1.D0 - rDensity(3) ) ) * d2rDensity_d2Volume(3) ) - ( ( 1.D0 / ( 1.D0 - &
&     rDensity(3) ) ) * drDensity_dVolume(3) * ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 2.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) &
&     * rDensity(3) ) ) - drDensity_dVolume(0) ) ) + ( DLOG( 1.D0 - rDensity(3) ) * ( ( ( 3.D0 * ( ( 2.D0 * rDensity(2) * &
&     drDensity_dVolume(2) * drDensity_dVolume(2) ) + ( rDensity(2) * rDensity(2) * d2rDensity_d2Volume(2) ) ) ) / ( rDensity(3) * &
&     rDensity(3) ) ) - ( ( 6.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) * drDensity_dVolume(3) ) / ( rDensity(3) * &
&     rDensity(3) * rDensity(3) ) ) - ( ( 2.D0 * ( ( 3.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) * &
&     drDensity_dVolume(3) ) + ( d2rDensity_d2Volume(3) * rDensity(2) * rDensity(2) * rDensity(2) ) ) ) / ( rDensity(3) * &
&     rDensity(3) * rDensity(3) ) ) + ( ( 6.D0 * rDensity(2) * rDensity(2) * rDensity(2) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) / ( rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) ) ) - d2rDensity_d2Volume(0) ) )
d2AuxiliaryHSFactor_d2Volume(2) = ( ( 3.D0 * ( ( d2rDensity_d2Volume(1) * rDensity(2) ) + ( 2.D0 * drDensity_dVolume(1) * &
&     drDensity_dVolume(2) ) + ( d2rDensity_d2Volume(2) * rDensity(1) ) ) ) / ( 1.D0 - rDensity(3) ) ) + ( ( 3.D0 * &
&     ( ( drDensity_dVolume(1) * rDensity(2) ) + ( drDensity_dVolume(2) * rDensity(1) ) ) * drDensity_dVolume(3) ) / ( ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 3.D0 * ( ( drDensity_dVolume(1) * drDensity_dVolume(3) * rDensity(2) ) + &
&     ( drDensity_dVolume(2) * drDensity_dVolume(3) * rDensity(1) ) + ( rDensity(1) * rDensity(2) * d2rDensity_d2Volume(3) ) ) ) / &
&     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 6.D0 * rDensity(1) * rDensity(2) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
d2AuxiliaryHSFactor_d2Volume(3) = - ( ( ( 2.D0 * rDensity(2) * drDensity_dVolume(2) * ( 3.D0 * drDensity_dVolume(2) * rDensity(3) &
&     * rDensity(3) - ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * drDensity_dVolume(2) ) * rDensity(3) + &
&     drDensity_dVolume(3) * rDensity(2) ) ) + ( rDensity(2) * rDensity(2) * ( ( 6.D0 * rDensity(3) * drDensity_dVolume(2) * &
&     drDensity_dVolume(3) ) + ( 3.D0 * rDensity(3) * rDensity(3) * d2rDensity_d2Volume(2) ) + ( d2rDensity_d2Volume(3) * &
&     rDensity(2) ) + ( drDensity_dVolume(3) * drDensity_dVolume(2) ) - ( drDensity_dVolume(3) * ( 3.D0 * drDensity_dVolume(3) * &
&     rDensity(2) + 3.D0 * drDensity_dVolume(2) ) ) - ( 3.D0 * rDensity(3) * ( ( d2rDensity_d2Volume(3) * rDensity(2) ) + &
&     ( drDensity_dVolume(3) * drDensity_dVolume(2) ) + d2rDensity_d2Volume(2) ) ) ) ) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) * rDensity(3) ) ) + ( ( ( rDensity(2) * rDensity(2) * ( 3.D0 * &
&     drDensity_dVolume(2) * rDensity(3) * rDensity(3) - ( 3.D0 * drDensity_dVolume(3) * rDensity(2) + 3.D0 * &
&     drDensity_dVolume(2) ) * rDensity(3) + drDensity_dVolume(3) * rDensity(2) ) ) * ( drDensity_dVolume(3) * ( 2.D0 - ( 5.D0 * &
&     rDensity(3) ) ) ) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     rDensity(3) * rDensity(3) * rDensity(3) ) )

! Second derivative of the Boublik-Mansoori hard-sphere mixture term with respect to the volume [mol² / m⁶]
d2HSBoublikMansoori_d2Volume = ( 12.D0 / cPi ) * dInverseNumberDensity_dVolume * SUM( dAuxiliaryHSFactor_dVolume ) + ( 6.D0 / &
&     ( cPi * rNumberDensity ) ) * SUM( d2AuxiliaryHSFactor_d2Volume )

! Second derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2ExcludedVolumeFreeEnergy_d2Volume = d2HSBoublikMansoori_d2Volume * AngleAverage * cUniversalGas * Temperature ! Proven units
IF( ReferenceBoublikLogical ) THEN
  d2ExcludedVolumeFreeEnergy_d2Volume = ( ( 2.D0 * ( 1.D0 + 3.D0 * cNonSphericity ) * rDensity(3) ) + ( ( 9.D0 * cNonSphericity * &
  &     cNonSphericity - 6.D0 * cNonSphericity - 5.D0 ) * rDensity(3) * rDensity(3) ) + ( ( 1.D0 - cNonSphericity * &
  &     cNonSphericity ) * ( 4.D0 - rDensity(3) ) * rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * mVolume * &
  &     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) )
  d2ExcludedVolumeFreeEnergy_d2Volume = d2ExcludedVolumeFreeEnergy_d2Volume * cUniversalGas * Temperature
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! Second derivative of the effective packing fraction with respect to the volume [mol² / m⁶]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  d2EffPackingFraction_d2Volume = - dEffPackingFraction_dVolume / mVolume + ( EffPackingFractionCoefficients(1) * rDensity(3) + &
  &     4.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) + 9.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     EffPackingFractionCoefficients(3) ) / ( mVolume * mVolume )
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2EffPackingFraction_d2Volume = - dEffPackingFraction_dVolume / mVolume + ( EffPackingFractionCoefficients(1) * rDensity(3) + &
  &     4.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) ) / ( mVolume * mVolume )
  d2EffPackingFractionSL_d2Volume = - dEffPackingFractionSL_dVolume / mVolume + ( EffPackingFractionCoefficientsSutherland(1) * &
  &     rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsSutherland(2) ) / ( mVolume * mVolume )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2EffPackingFraction_d2Volume = - dEffPackingFraction_dVolume / mVolume + ( EffPackingFractionCoefficients(1) * rDensity(3) + &
  &     4.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) ) / ( mVolume * mVolume )
  d2EffPackingFractionYK_d2Volume = - dEffPackingFractionYK_dVolume / mVolume + ( EffPackingFractionCoefficientsYukawa(1) * &
  &     rDensity(3) + 4.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsYukawa(2) ) / ( mVolume * mVolume )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2EffPackingFraction_d2Volume = ( ( 2.D0 * EffPackingFractionCoefficients(1) * rDensity(3) ) + ( 6.D0 * &
  &     EffPackingFractionCoefficients(2) * rDensity(3) * rDensity(3) ) + ( 12.D0 * EffPackingFractionCoefficients(3) * &
  &     rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * mVolume )
  d2EffPackingFractionCSW_d2Volume = ( ( 2.D0 * EffPackingFractionCoefficientsCSW(1) * rDensity(3) ) + ( 6.D0 * &
  &     EffPackingFractionCoefficientsCSW(2) * rDensity(3) * rDensity(3) ) + ( 12.D0 * EffPackingFractionCoefficientsCSW(3) * &
  &     rDensity(3) * rDensity(3) * rDensity(3) ) ) / ( mVolume * mVolume )
END IF

! First derivative with respect to the volume of the cross derivative of the effective packing fraction with respect to the volume and the density [mol² . (Å³ / m⁶)]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  d2EffPFraction_d2Volume_dDensity = dxEffPFraction_dVolume_dDensity / mVolume - dEffPackingFraction_dDensity / ( mVolume * &
  &     mVolume ) + ( 1.D0 / ( rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1) * rDensity(3) ) + &
  &     ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) ) + ( 27.D0 * rDensity(3) * rDensity(3) * &
  &     rDensity(3) * EffPackingFractionCoefficients(3) ) )
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2EffPFraction_d2Volume_dDensity = dxEffPFraction_dVolume_dDensity / mVolume - dEffPackingFraction_dDensity / ( mVolume * &
  &     mVolume ) + ( 1.D0 / ( rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1) * rDensity(3) ) + &
  &     ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) ) )
  d2EffPFractionSL_d2Volume_dDensity = dxEffPFractionSL_dVolume_dDensity / mVolume - dEffPackingFractionSutherland_dDensity / &
  &     ( mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsSutherland(1) &
  &     * rDensity(3) ) + ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsSutherland(2) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2EffPFraction_d2Volume_dDensity = dxEffPFraction_dVolume_dDensity / mVolume - dEffPackingFraction_dDensity / ( mVolume * &
  &     mVolume ) + ( 1.D0 / ( rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficients(1) * rDensity(3) ) + &
  &     ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(2) ) )
  d2EffPFractionYK_d2Volume_dDensity = dxEffPFractionYK_dVolume_dDensity / mVolume - dEffPackingFractionYukawa_dDensity / ( &
  &     mVolume * mVolume ) + ( 1.D0 / ( rNumberDensity * mVolume * mVolume ) ) * ( (EffPackingFractionCoefficientsYukawa(1) * &
  &     rDensity(3) ) + ( 8.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsYukawa(2) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2EffPFraction_d2Volume_dDensity = ( 4.D0 * EffPackingFractionCoefficients(2) * rDensity(3) * rDensity(3) / ( rNumberDensity * &
  &     mVolume * mVolume ) ) + ( 18.D0 * EffPackingFractionCoefficients(3) * rDensity(3) * rDensity(3) * rDensity(3) / &
  &     ( rNumberDensity * mVolume * mVolume ) )
  d2EffPFractionCSW_d2Volume_dDensity = ( 4.D0 * EffPackingFractionCoefficientsCSW(2) * rDensity(3) * rDensity(3) / &
  &     ( rNumberDensity * mVolume * mVolume ) ) + ( 18.D0 * EffPackingFractionCoefficientsCSW(3) * rDensity(3) * rDensity(3) * &
  &     rDensity(3) / ( rNumberDensity * mVolume * mVolume ) )
END IF

! Second derivative of the contact radial distribution function with respect to the volume [mol² / m⁶]
d2cRDFunction_d2Volume = ( ( d2rDensity_d2Volume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) * ( 1.D0 + ( 3.D0 * &
&     AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + 2.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) + ( 2.D0 * &
&     AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
&     rDensity(3) / ( 1.D0 - rDensity(3) ) ) + ( 3.D0 * rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) ) ) ) ) ) + ( ( 2.D0 * drDensity_dVolume(3) * drDensity_dVolume(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) * ( 1.D0 + ( 3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship * &
&     ( 2.D0 + ( 3.D0 * rDensity(3) / ( 1.D0 - rDensity(3) ) ) ) ) + ( 2.D0 * cAverageDiameterRelationship * &
&     cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship * ( 1.D0 + ( 6.D0 * rDensity(3) / ( 1.D0 - &
&     rDensity(3) ) ) + ( 6.D0 * rDensity(3) * rDensity(3) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) ) ) ) )

! Second derivative of the radial distribution function for an effective packing fraction with respect to the volume [mol² / m⁶]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  d2cEffectiveRDF_d2Volume = ( ( d2EffPackingFraction_d2Volume / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) &
  &     ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + 2.D0 * EffPackingFraction / ( 1.D0 &
  &     - EffPackingFraction ) ) ) + ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * &
  &     cAverageDiameterRelationship * ( ( 2.D0 * EffPackingFraction / ( 1.D0 - EffPackingFraction ) ) + ( 3.D0 * &
  &     EffPackingFraction * EffPackingFraction / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) ) ) ) ) + &
  &     ( ( ( 2.D0 * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) * ( 1.D0 + ( ( 3.D0 * AuxDiameterRelationship * &
  &     cAverageDiameterRelationship ) * ( 2.D0 + ( 3.D0 * EffPackingFraction / ( 1.D0 - EffPackingFraction ) ) ) ) + ( ( 2.D0 * &
  &     AuxDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship ) * ( 1.D0 &
  &     + ( 6.D0 * EffPackingFraction / ( 1.D0 - EffPackingFraction ) ) + ( 6.D0 * EffPackingFraction * EffPackingFraction / &
  &     ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) ) ) ) )
ELSE
  d2cEffectiveRDF_d2Volume = d2EffPackingFraction_d2Volume / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) + &
  &     ( ( dEffPackingFraction_dVolume * dEffPackingFraction_dVolume * ( 2.D0 + 3.D0 * gAux1 ) ) + ( 3.D0 * gAux1 * ( 1.D0 + &
  &     EffPackingFraction ) * d2EffPackingFraction_d2Volume ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) ) + ( ( ( 4.D0 * gAux2 + 9.D0 * gAux1 ) * ( 1.D0 + EffPackingFraction ) * &
  &     dEffPackingFraction_dVolume * dEffPackingFraction_dVolume ) + ( 2.D0 * gAux2 * EffPackingFraction * ( 2.D0 + &
  &     EffPackingFraction ) * d2EffPackingFraction_d2Volume ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) + ( 8.D0 * gAux2 * EffPackingFraction * ( 2.D0 + &
  &     EffPackingFraction ) * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume ) / ( ( 1.D0 - EffPackingFraction ) * &
  &     ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2cEffectiveRDFSutherland_d2Volume = ( ( d2EffPackingFractionSL_d2Volume / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + &
  &     2.D0 * EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( 2.D0 * AuxDiameterRelationship * &
  &     cAverageDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
  &     EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) + ( 3.D0 * EffPackingFractionSutherland * &
  &     EffPackingFractionSutherland / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) ) ) &
  &     ) ) + ( ( ( 2.D0 * dEffPackingFractionSL_dVolume * dEffPackingFractionSL_dVolume ) / ( ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) * ( &
  &     1.D0 + ( ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship ) * ( 2.D0 + ( 3.D0 * &
  &     EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) ) ) + ( ( 2.D0 * AuxDiameterRelationship * &
  &     AuxDiameterRelationship * cAverageDiameterRelationship * cAverageDiameterRelationship ) * ( 1.D0 + ( 6.D0 * &
  &     EffPackingFractionSutherland / ( 1.D0 - EffPackingFractionSutherland ) ) + ( 6.D0 * EffPackingFractionSutherland * &
  &     EffPackingFractionSutherland / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) ) ) &
  &     ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2cEffectiveRDFYukawa_d2Volume = ( ( d2EffPackingFractionYK_d2Volume / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) * ( 1.D0 + ( 3.D0 * AuxDiameterRelationship * cAverageDiameterRelationship * ( 1.D0 + &
  &     2.D0 * EffPackingFractionYukawa / ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( 2.D0 * AuxDiameterRelationship * &
  &     cAverageDiameterRelationship * AuxDiameterRelationship * cAverageDiameterRelationship * ( ( 2.D0 * &
  &     EffPackingFractionYukawa / ( 1.D0 - EffPackingFractionYukawa ) ) + ( 3.D0 * EffPackingFractionYukawa * &
  &     EffPackingFractionYukawa / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) ) ) ) ) + &
  &     ( ( ( 2.D0 * dEffPackingFractionYK_dVolume * dEffPackingFractionYK_dVolume ) / ( ( 1.D0 - EffPackingFractionYukawa ) * &
  &     ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) * ( 1.D0 + ( ( 3.D0 * &
  &     AuxDiameterRelationship * cAverageDiameterRelationship ) * ( 2.D0 + ( 3.D0 * EffPackingFractionYukawa / ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) ) + ( ( 2.D0 * AuxDiameterRelationship * AuxDiameterRelationship * &
  &     cAverageDiameterRelationship * cAverageDiameterRelationship ) * ( 1.D0 + ( 6.D0 * EffPackingFractionYukawa / ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) + ( 6.D0 * EffPackingFractionYukawa * EffPackingFractionYukawa / ( ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) ) ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2cEffectiveRDFCSW_d2Volume = d2EffPackingFractionCSW_d2Volume / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) + ( ( dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dVolume * ( 2.D0 + 3.D0 * &
  &     gAux1 ) ) + ( 3.D0 * gAux1 * ( 1.D0 + EffPackingFractionCSW ) * d2EffPackingFractionCSW_d2Volume ) ) / ( ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) + ( ( ( 4.D0 * gAux2 + &
  &     9.D0 * gAux1 ) * ( 1.D0 + EffPackingFractionCSW ) * dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dVolume ) + &
  &     ( 2.D0 * gAux2 * EffPackingFractionCSW * ( 2.D0 + EffPackingFractionCSW ) * d2EffPackingFractionCSW_d2Volume ) ) / ( ( &
  &     1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) + ( 8.D0 * gAux2 * EffPackingFractionCSW * ( 2.D0 + EffPackingFractionCSW ) * &
  &     dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dVolume ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) )
END IF

! Second derivative of the mean-attractive energy with respect to the volume [K . mol² / m⁶]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  d2cMeanAttractiveEnergy_d2Volume = - 2.D0 * cPi * ( ijPotentialRangeCubic(cComponent,cComponent) - 1.D0 ) * &
  &     ( ijaDiameterSphereCubic(cComponent,cComponent) / 3.D0 ) * ( ( d2NumberDensity_d2Volume * &
  &     cEffectiveRadialDistributionFunction ) + ( 2.D0 * dcEffectiveRDF_dVolume * dNumberDensity_dVolume ) + ( rNumberDensity * &
  &     d2cEffectiveRDF_d2Volume ) ) * ijaWellDepth(cComponent,cComponent)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2cMeanAttractiveEnergy_d2Volume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) * &
  &     ( ( d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunction ) + ( 2.D0 * dcEffectiveRDF_dVolume * &
  &     dNumberDensity_dVolume ) + ( rNumberDensity * d2cEffectiveRDF_d2Volume ) )
  d2cMeanAttractiveEnergySL_d2Volume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) * &
  &     ( ( d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunctionSL ) + ( 2.D0 * dcEffectiveRDFSutherland_dVolume * &
  &     dNumberDensity_dVolume ) + ( rNumberDensity * d2cEffectiveRDFSutherland_d2Volume ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2cMeanAttractiveEnergy_d2Volume = - 2.D0 * cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) + ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) ) ) * ( ( d2NumberDensity_d2Volume * &
  &     cEffectiveRadialDistributionFunction ) + ( 2.D0 * dcEffectiveRDF_dVolume * dNumberDensity_dVolume ) + ( rNumberDensity * &
  &     d2cEffectiveRDF_d2Volume ) )
  d2cMeanAttractiveEnergyYK_d2Volume = - cPi * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) * ( ( &
  &     d2NumberDensity_d2Volume * cEffectiveRadialDistributionFunctionYK ) + ( 2.D0 * dcEffectiveRDFYukawa_dVolume * &
  &     dNumberDensity_dVolume ) + ( rNumberDensity * d2cEffectiveRDFYukawa_d2Volume ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2cMeanAttractiveEnergy_d2Volume = - ijaWellDepth(cComponent,cComponent) * ( fMolecularVolume * ( 1.D0 + 3.D0 * fNonSphericity ) &
  &     - cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) ) * ( ( cEffectiveRadialDistributionFunction * &
  &     d2NumberDensity_d2Volume ) + ( 2.D0 * dNumberDensity_dVolume * dcEffectiveRDF_dVolume ) + ( d2cEffectiveRDF_d2Volume * &
  &     rNumberDensity ) )
  d2cMeanAttractiveEnergyCSW_d2Volume = - ijaWellDepth(cComponent,cComponent) * ( fMolecularVolume * ( 1.D0 + 3.D0 * &
  &     fNonSphericity ) - cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) ) * ( ( cEffectiveRadialDistributionFunctionCSW * &
  &     d2NumberDensity_d2Volume ) + ( 2.D0 * dNumberDensity_dVolume * dcEffectiveRDFCSW_dVolume ) + ( rNumberDensity * &
  &     d2cEffectiveRDFCSW_d2Volume ) )
END IF

! Second derivative of the first-order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MeanAttractiveFreeEnergy_d2Volume = d2cMeanAttractiveEnergy_d2Volume * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Auxiliary factors of a second derivative of the hard-sphere isothermal compressibility with respect to the volume
d2AuxHSIC_d2Volume(1) = ( ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) ) - ( 4.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     drDensity_dVolume(3) * rDensity(0) ) ) ! [mol / (Å³ . m³)]
d2AuxHSIC_d2Volume(2) = ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * &
&     ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) ! [1 / Å³]
d2AuxHSIC_d2Volume(3) = ( d2rDensity_d2Volume(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) ) ) - ( 8.D0 * drDensity_dVolume(0) * drDensity_dVolume(3) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( 12.D0 * rDensity(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) * &
&     drDensity_dVolume(3) * drDensity_dVolume(3) ) - ( 4.D0 * rDensity(0) * d2rDensity_d2Volume(3) * ( ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) ! [mol² / (Å³ . m⁶)]
d2AuxHSIC_d2Volume(4) = ( drDensity_dVolume(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * rDensity(0) * drDensity_dVolume(3) ) + ( 6.D0 * ( ( drDensity_dVolume(1) * rDensity(2) * ( 1.D0 - &
&     rDensity(3) ) ) + ( drDensity_dVolume(2) * rDensity(1) * ( 1.D0 - rDensity(3) ) ) - ( drDensity_dVolume(3) * rDensity(1) * &
&     rDensity(2) ) ) ) + ( 27.D0 * rDensity(2) * rDensity(2) * drDensity_dVolume(2) ) ! [mol / (Å³ . m³)]
d2AuxHSIC_d2Volume(5) = ( d2AuxHSIC_d2Volume(1) / ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) - ( ( ( 2.D0 * rDensity(0) * &
&     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) / &
&     ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * d2AuxHSIC_d2Volume(4) ) ! [mol . (Å³ / m³)]
d2AuxHSIC_d2Volume(6) = - ( ( rDensity(0) * ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) ) / ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * ( ( d2rDensity_d2Volume(0) * ( ( 1.D0 - &
&     rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( 4.D0 * drDensity_dVolume(0) * drDensity_dVolume(3) * ( 1.D0 - rDensity(3) ) ) &
&     - ( 2.D0 * rDensity(0) * ( 1.D0 - rDensity(3) ) * d2rDensity_d2Volume(3) ) + ( 2.D0 * rDensity(0) * drDensity_dVolume(3) * &
&     drDensity_dVolume(3) ) + ( 6.D0 * ( ( d2rDensity_d2Volume(1) * rDensity(2) * ( 1.D0 - rDensity(3) ) ) + ( 2.D0 * ( 1.D0 - &
&     rDensity(3) ) * drDensity_dVolume(1) * drDensity_dVolume(2) ) - ( 2.D0 * drDensity_dVolume(1) * drDensity_dVolume(3) * &
&     rDensity(2) ) + ( rDensity(1) * d2rDensity_d2Volume(2) * ( 1.D0 - rDensity(3) ) ) - ( 2.D0 * rDensity(1) * &
&     drDensity_dVolume(2) * drDensity_dVolume(3) ) - ( rDensity(1) * rDensity(2) * d2rDensity_d2Volume(3) ) ) ) + ( 54.D0 * &
&     rDensity(2) * drDensity_dVolume(2) * drDensity_dVolume(2) ) + ( 27.D0 * rDensity(2) * rDensity(2) * &
&     d2rDensity_d2Volume(2) ) ) ! [mol² / m⁶]

! Second derivative of the hard-sphere isothermal compressibility with respect to the volume [mol² / m⁶]
d2HSIsothermalCompressibility_d2Volume = ( d2AuxHSIC_d2Volume(3) / d2AuxHSIC_d2Volume(2) ) - ( ( d2AuxHSIC_d2Volume(1) / &
&     ( d2AuxHSIC_d2Volume(2) * d2AuxHSIC_d2Volume(2) ) ) * d2AuxHSIC_d2Volume(4) ) - ( d2AuxHSIC_d2Volume(5) * &
&     d2AuxHSIC_d2Volume(4) ) + d2AuxHSIC_d2Volume(6)

! Second derivative of the hard convex-body isothermal compressibility with respect to the volume [mol² / m⁶]
IF( PYHCBCorrectionLogical ) THEN
  d2AuxHCBIC_d2Volume(1) = 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 )
  d2AuxHCBIC_d2Volume(2) = ( 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + 2.D0 * rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) ) / mVolume
  d2AuxHCBIC_d2Volume(3) = - ( 4.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + 6.D0 * rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) ) / ( mVolume * mVolume )
  d2HCBIsothermalCompressibility_d2Volume = ( ( ( 12.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) * &
  &     rDensity(3) ) - ( 8.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) ) ) / &
  &     ( mVolume * mVolume * d2AuxHCBIC_d2Volume(1) ) ) + ( ( ( d2AuxHCBIC_d2Volume(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 8.D0 * d2AuxHCBIC_d2Volume(2) * ( rDensity(3) / &
  &     mVolume ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) / ( d2AuxHCBIC_d2Volume(1) * &
  &     d2AuxHCBIC_d2Volume(1) ) ) + ( ( 2.D0 * d2AuxHCBIC_d2Volume(2) * d2AuxHCBIC_d2Volume(2) * ( 1.D0 - rDensity(3) ) * &
  &     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) / ( d2AuxHCBIC_d2Volume(1) * &
  &     d2AuxHCBIC_d2Volume(1) * d2AuxHCBIC_d2Volume(1) ) )
ELSE
  d2AuxHCBIC_d2Volume(1) = 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 4.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     cNonSphericity * cNonSphericity + rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * cNonSphericity
  d2AuxHCBIC_d2Volume(2) = ( 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + 2.D0 * rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 12.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     cNonSphericity * cNonSphericity + 4.D0 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity ) / mVolume
  d2AuxHCBIC_d2Volume(3) = - ( 4.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + 6.D0 * rDensity(3) * rDensity(3) * ( 9.D0 * &
  &     cNonSphericity * cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 48.D0 * rDensity(3) * rDensity(3) * rDensity(3) * &
  &     cNonSphericity * cNonSphericity + 20.D0 * rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity ) / ( mVolume * mVolume )
  d2HCBIsothermalCompressibility_d2Volume = ( ( ( 12.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) * &
  &     rDensity(3) ) - ( 8.D0 * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * rDensity(3) ) ) / &
  &     ( mVolume * mVolume * d2AuxHCBIC_d2Volume(1) ) ) + ( ( ( d2AuxHCBIC_d2Volume(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
  &     rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 8.D0 * d2AuxHCBIC_d2Volume(2) * ( rDensity(3) / &
  &     mVolume ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) / ( d2AuxHCBIC_d2Volume(1) * &
  &     d2AuxHCBIC_d2Volume(1) ) ) + ( ( 2.D0 * d2AuxHCBIC_d2Volume(2) * d2AuxHCBIC_d2Volume(2) * ( 1.D0 - rDensity(3) ) * &
  &     ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) / ( d2AuxHCBIC_d2Volume(1) * &
  &     d2AuxHCBIC_d2Volume(1) * d2AuxHCBIC_d2Volume(1) ) )
END IF

! Second derivative of the Zhang's factor with respect to the volume [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    d2ZhangFactor_d2Volume = 4.D0 * ZhangCorrection * ( ( drDensity_dVolume(3) * drDensity_dVolume(3) ) + ( rDensity(3) * &
    &     d2rDensity_d2Volume(3) ) )
  ELSE
    d2ZhangFactor_d2Volume = 0.75D0 * ZhangCorrection * ijSecondVirialCoefficient(cComponent,cComponent) * &
    &     ijSecondVirialCoefficient(cComponent,cComponent) * rNumberDensity * rNumberDensity / ( mVolume * mVolume )
  END IF
END IF

! First derivative with respect to the volume of an auxiliary factor of the cross derivative of the effective radial distribution function with respect to the volume and the density [mol² . (Å³ / m⁶)]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  d2AuxEffRDF_d2Volume_dDensity(1) = ( ( d2EffPFraction_d2Volume_dDensity * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) + ( 2.D0 * ( 1.D0 - EffPackingFraction ) * dxEffPFraction_dVolume_dDensity * &
  &     dEffPackingFraction_dVolume ) + ( 6.D0 * dEffPackingFraction_dDensity * dEffPackingFraction_dVolume * &
  &     dEffPackingFraction_dVolume ) + ( 2.D0 * ( 1.D0 - EffPackingFraction ) * ( ( dEffPackingFraction_dVolume * &
  &     dxEffPFraction_dVolume_dDensity ) + ( dEffPackingFraction_dDensity * d2EffPackingFraction_d2Volume ) ) ) ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) )
  d2AuxEffRDF_d2Volume_dDensity(2) = 3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship * ( ( ( ( ( 1.D0 - &
  &     EffPackingFraction * EffPackingFraction ) * d2EffPFraction_d2Volume_dDensity ) - ( 2.D0 * dxEffPFraction_dVolume_dDensity &
  &     * EffPackingFraction * dEffPackingFraction_dVolume ) + ( 2.D0 * dxEffPFraction_dVolume_dDensity * &
  &     dEffPackingFraction_dVolume * ( 2.D0 + EffPackingFraction ) ) + ( 2.D0 * d2EffPackingFraction_d2Volume * &
  &     dEffPackingFraction_dDensity * ( 2.D0 + EffPackingFraction ) ) + ( 2.D0 * dEffPackingFraction_dDensity * &
  &     dEffPackingFraction_dVolume * dEffPackingFraction_dVolume ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + ( ( ( 4.D0 * &
  &     dxEffPFraction_dVolume_dDensity * dEffPackingFraction_dVolume * ( 1.D0 - EffPackingFraction * EffPackingFraction ) ) + &
  &     ( 8.D0 * dEffPackingFraction_dDensity * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume * ( 2.D0 + &
  &     EffPackingFraction ) ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) )
  d2AuxEffRDF_d2Volume_dDensity(3) = 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * AuxDiameterRelationship &
  &     * AuxDiameterRelationship * ( ( ( ( 2.D0 * dxEffPFraction_dVolume_dDensity * dEffPackingFraction_dVolume * ( 1.D0 + 4.D0 * &
  &     EffPackingFraction + EffPackingFraction * EffPackingFraction ) ) + ( 2.D0 * dEffPackingFraction_dDensity * &
  &     d2EffPackingFraction_d2Volume * ( 1.D0 + 4.D0 * EffPackingFraction + EffPackingFraction * EffPackingFraction ) ) + ( 4.D0 &
  &     * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * ( 2.D0 + EffPackingFraction ) &
  &     ) + ( d2EffPFraction_d2Volume_dDensity * EffPackingFraction * ( 2.D0 - EffPackingFraction - EffPackingFraction * &
  &     EffPackingFraction ) ) + ( dxEffPFraction_dVolume_dDensity * dEffPackingFraction_dVolume * ( 2.D0 - EffPackingFraction - &
  &     EffPackingFraction * EffPackingFraction ) ) - ( dxEffPFraction_dVolume_dDensity * EffPackingFraction * &
  &     dEffPackingFraction_dVolume * ( 1.D0 + 2.D0 * EffPackingFraction ) ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + &
  &     ( ( 5.D0 * dEffPackingFraction_dVolume * ( ( 2.D0 * dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * ( 1.D0 + &
  &     4.D0 * EffPackingFraction + EffPackingFraction * EffPackingFraction ) ) + ( dxEffPFraction_dVolume_dDensity * &
  &     EffPackingFraction * ( 2.D0 - EffPackingFraction - EffPackingFraction * EffPackingFraction ) ) ) ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * &
  &     ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) )
ELSE
  d2AuxEffRDF_d2Volume_dDensity(1) = d2EffPFraction_d2Volume_dDensity / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) + ( ( 3.D0 * gAux1 + 2.D0 ) * ( 2.D0 * dxEffPFraction_dVolume_dDensity * &
  &     dEffPackingFraction_dVolume + dEffPackingFraction_dDensity * d2EffPackingFraction_d2Volume ) + 3.D0 * gAux1 * &
  &     d2EffPFraction_d2Volume_dDensity * ( 1.D0 + EffPackingFraction ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) )
  d2AuxEffRDF_d2Volume_dDensity(2) = ( 2.D0 * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume * &
  &     dEffPackingFraction_dDensity * ( 9.D0 * gAux1 + 2.D0 * gAux2 + 3.D0 ) + dEffPackingFraction_dDensity * &
  &     d2EffPackingFraction_d2Volume * ( 1.D0 + EffPackingFraction ) * ( 9.D0 * gAux1 + 4.D0 * gAux2 ) + 2.D0 * &
  &     dxEffPFraction_dVolume_dDensity * dEffPackingFraction_dVolume * ( 1.D0 + EffPackingFraction ) * ( 9.D0 * gAux1 + 4.D0 * &
  &     gAux2 ) + 2.D0 * d2EffPFraction_d2Volume_dDensity * gAux2 * EffPackingFraction * ( 2.D0 + EffPackingFraction ) ) / &
  &     ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) )
  d2AuxEffRDF_d2Volume_dDensity(3) = ( 16.D0 * dxEffPFraction_dVolume_dDensity * dEffPackingFraction_dVolume * EffPackingFraction &
  &     * ( 2.D0 + EffPackingFraction ) * gAux2 + 4.D0 * dEffPackingFraction_dVolume * dEffPackingFraction_dVolume * &
  &     dEffPackingFraction_dDensity * ( 1.D0 + EffPackingFraction ) * ( 9.D0 * gAux1 + 8.D0 * gAux2 ) + 8.D0 * &
  &     dEffPackingFraction_dDensity * d2EffPackingFraction_d2Volume * EffPackingFraction * ( 2.D0 + EffPackingFraction ) * &
  &     gAux2 ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) + ( 40.D0 * gAux2 * dEffPackingFraction_dVolume * &
  &     dEffPackingFraction_dVolume * dEffPackingFraction_dDensity * EffPackingFraction * ( 2.D0 + EffPackingFraction ) ) / &
  &     ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2AuxEffRDFSutherland_d2Volume_dDensity(1) = ( ( d2EffPFractionSL_d2Volume_dDensity * ( 1.D0 - EffPackingFractionSutherland ) * &
  &     ( 1.D0 - EffPackingFractionSutherland ) ) + ( 2.D0 * ( 1.D0 - EffPackingFractionSutherland ) * &
  &     dxEffPFractionSL_dVolume_dDensity * dEffPackingFractionSL_dVolume ) + ( 6.D0 * dEffPackingFractionSutherland_dDensity * &
  &     dEffPackingFractionSL_dVolume * dEffPackingFractionSL_dVolume ) + ( 2.D0 * ( 1.D0 - EffPackingFractionSutherland ) * ( ( &
  &     dEffPackingFractionSL_dVolume * dxEffPFractionSL_dVolume_dDensity ) + ( dEffPackingFractionSutherland_dDensity * &
  &     d2EffPackingFractionSL_d2Volume ) ) ) ) / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) )
  d2AuxEffRDFSutherland_d2Volume_dDensity(2) = 3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship * ( ( ( ( ( 1.D0 - &
  &     EffPackingFractionSutherland * EffPackingFractionSutherland ) * d2EffPFractionSL_d2Volume_dDensity ) - ( 2.D0 * &
  &     dxEffPFractionSL_dVolume_dDensity * EffPackingFractionSutherland * dEffPackingFractionSL_dVolume ) + ( 2.D0 * &
  &     dxEffPFractionSL_dVolume_dDensity * dEffPackingFractionSL_dVolume * ( 2.D0 + EffPackingFractionSutherland ) ) + ( 2.D0 * &
  &     d2EffPackingFractionSL_d2Volume * dEffPackingFractionSutherland_dDensity * ( 2.D0 + EffPackingFractionSutherland ) ) + ( &
  &     2.D0 * dEffPackingFractionSutherland_dDensity * dEffPackingFractionSL_dVolume * dEffPackingFractionSL_dVolume ) ) / ( ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) &
  &     * ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( ( ( 4.D0 * dxEffPFractionSL_dVolume_dDensity * &
  &     dEffPackingFractionSL_dVolume * ( 1.D0 - EffPackingFractionSutherland * EffPackingFractionSutherland ) ) + ( 8.D0 * &
  &     dEffPackingFractionSutherland_dDensity * dEffPackingFractionSL_dVolume * dEffPackingFractionSL_dVolume * ( 2.D0 + &
  &     EffPackingFractionSutherland ) ) ) / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * &
  &     ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland &
  &     ) ) ) )
  d2AuxEffRDFSutherland_d2Volume_dDensity(3) = 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * &
  &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 * dxEffPFractionSL_dVolume_dDensity * &
  &     dEffPackingFractionSL_dVolume * ( 1.D0 + 4.D0 * EffPackingFractionSutherland + EffPackingFractionSutherland * &
  &     EffPackingFractionSutherland ) ) + ( 2.D0 * dEffPackingFractionSutherland_dDensity * d2EffPackingFractionSL_d2Volume * ( &
  &     1.D0 + 4.D0 * EffPackingFractionSutherland + EffPackingFractionSutherland * EffPackingFractionSutherland ) ) + ( 4.D0 * &
  &     dEffPackingFractionSL_dVolume * dEffPackingFractionSL_dVolume * dEffPackingFractionSutherland_dDensity * ( 2.D0 + &
  &     EffPackingFractionSutherland ) ) + ( d2EffPFractionSL_d2Volume_dDensity * EffPackingFractionSutherland * ( 2.D0 - &
  &     EffPackingFractionSutherland - EffPackingFractionSutherland * EffPackingFractionSutherland ) ) + ( &
  &     dxEffPFractionSL_dVolume_dDensity * dEffPackingFractionSL_dVolume * ( 2.D0 - EffPackingFractionSutherland - &
  &     EffPackingFractionSutherland * EffPackingFractionSutherland ) ) - ( dxEffPFractionSL_dVolume_dDensity * &
  &     EffPackingFractionSutherland * dEffPackingFractionSL_dVolume * ( 1.D0 + 2.D0 * EffPackingFractionSutherland ) ) ) / ( ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) &
  &     * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( ( 5.D0 * &
  &     dEffPackingFractionSL_dVolume * ( ( 2.D0 * dEffPackingFractionSL_dVolume * dEffPackingFractionSutherland_dDensity * ( 1.D0 &
  &     + 4.D0 * EffPackingFractionSutherland + EffPackingFractionSutherland * EffPackingFractionSutherland ) ) + &
  &     ( dxEffPFractionSL_dVolume_dDensity * EffPackingFractionSutherland * ( 2.D0 - EffPackingFractionSutherland - &
  &     EffPackingFractionSutherland * EffPackingFractionSutherland ) ) ) ) / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( &
  &     1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2AuxEffRDFYukawa_d2Volume_dDensity(1) = ( ( d2EffPFractionYK_d2Volume_dDensity * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) + ( 2.D0 * ( 1.D0 - EffPackingFractionYukawa ) * dxEffPFractionYK_dVolume_dDensity * &
  &     dEffPackingFractionYK_dVolume ) + ( 6.D0 * dEffPackingFractionYukawa_dDensity * dEffPackingFractionYK_dVolume * &
  &     dEffPackingFractionYK_dVolume ) + ( 2.D0 * ( 1.D0 - EffPackingFractionYukawa ) * ( ( dEffPackingFractionYK_dVolume * &
  &     dxEffPFractionYK_dVolume_dDensity ) + ( dEffPackingFractionYukawa_dDensity * d2EffPackingFractionYK_d2Volume ) ) ) ) / ( ( &
  &     1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) )
  d2AuxEffRDFYukawa_d2Volume_dDensity(2) = 3.D0 * cAverageDiameterRelationship * AuxDiameterRelationship * ( ( ( ( ( 1.D0 - &
  &     EffPackingFractionYukawa * EffPackingFractionYukawa ) * d2EffPFractionYK_d2Volume_dDensity ) - ( 2.D0 * &
  &     dxEffPFractionYK_dVolume_dDensity * EffPackingFractionYukawa * dEffPackingFractionYK_dVolume ) + ( 2.D0 * &
  &     dxEffPFractionYK_dVolume_dDensity * dEffPackingFractionYK_dVolume * ( 2.D0 + EffPackingFractionYukawa ) ) + ( 2.D0 * &
  &     d2EffPackingFractionYK_d2Volume * dEffPackingFractionYukawa_dDensity * ( 2.D0 + EffPackingFractionYukawa ) ) + ( 2.D0 * &
  &     dEffPackingFractionYukawa_dDensity * dEffPackingFractionYK_dVolume * dEffPackingFractionYK_dVolume ) ) / ( ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) ) ) + ( ( ( 4.D0 * dxEffPFractionYK_dVolume_dDensity * dEffPackingFractionYK_dVolume * ( 1.D0 - &
  &     EffPackingFractionYukawa * EffPackingFractionYukawa ) ) + ( 8.D0 * dEffPackingFractionYukawa_dDensity * &
  &     dEffPackingFractionYK_dVolume * dEffPackingFractionYK_dVolume * ( 2.D0 + EffPackingFractionYukawa ) ) ) / ( ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) )
  d2AuxEffRDFYukawa_d2Volume_dDensity(3) = 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * &
  &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( ( ( 2.D0 * dxEffPFractionYK_dVolume_dDensity * &
  &     dEffPackingFractionYK_dVolume * ( 1.D0 + 4.D0 * EffPackingFractionYukawa + EffPackingFractionYukawa * &
  &     EffPackingFractionYukawa ) ) + ( 2.D0 * dEffPackingFractionYukawa_dDensity * d2EffPackingFractionYK_d2Volume * ( 1.D0 + &
  &     4.D0 * EffPackingFractionYukawa + EffPackingFractionYukawa * EffPackingFractionYukawa ) ) + ( 4.D0 * &
  &     dEffPackingFractionYK_dVolume * dEffPackingFractionYK_dVolume * dEffPackingFractionYukawa_dDensity * ( 2.D0 + &
  &     EffPackingFractionYukawa ) ) + ( d2EffPFractionYK_d2Volume_dDensity * EffPackingFractionYukawa * ( 2.D0 - &
  &     EffPackingFractionYukawa - EffPackingFractionYukawa * EffPackingFractionYukawa ) ) + ( dxEffPFractionYK_dVolume_dDensity * &
  &     dEffPackingFractionYK_dVolume * ( 2.D0 - EffPackingFractionYukawa - EffPackingFractionYukawa * EffPackingFractionYukawa ) &
  &     ) - ( dxEffPFractionYK_dVolume_dDensity * EffPackingFractionYukawa * dEffPackingFractionYK_dVolume * ( 1.D0 + 2.D0 * &
  &     EffPackingFractionYukawa ) ) ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( ( 5.D0 * &
  &     dEffPackingFractionYK_dVolume * ( ( 2.D0 * dEffPackingFractionYK_dVolume * dEffPackingFractionYukawa_dDensity * ( 1.D0 + &
  &     4.D0 * EffPackingFractionYukawa + EffPackingFractionYukawa * EffPackingFractionYukawa ) ) + &
  &     ( dxEffPFractionYK_dVolume_dDensity * EffPackingFractionYukawa * ( 2.D0 - EffPackingFractionYukawa - &
  &     EffPackingFractionYukawa * EffPackingFractionYukawa ) ) ) ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2AuxEffRDFCSW_d2Volume_dDensity(1) = d2EffPFractionCSW_d2Volume_dDensity / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) + ( ( 3.D0 * gAux1 + 2.D0 ) * ( 2.D0 * dxEffPFractionCSW_dVolume_dDensity * &
  &     dEffPackingFractionCSW_dVolume + dEffPackingFractionCSW_dDensity * d2EffPackingFractionCSW_d2Volume ) + 3.D0 * gAux1 * &
  &     d2EffPFractionCSW_d2Volume_dDensity * ( 1.D0 + EffPackingFractionCSW ) ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) )
  d2AuxEffRDFCSW_d2Volume_dDensity(2) = ( 2.D0 * dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dVolume * &
  &     dEffPackingFractionCSW_dDensity * ( 9.D0 * gAux1 + 2.D0 * gAux2 + 3.D0 ) + dEffPackingFractionCSW_dDensity * &
  &     d2EffPackingFractionCSW_d2Volume * ( 1.D0 + EffPackingFractionCSW ) * ( 9.D0 * gAux1 + 4.D0 * gAux2 ) + 2.D0 * &
  &     dxEffPFractionCSW_dVolume_dDensity * dEffPackingFractionCSW_dVolume * ( 1.D0 + EffPackingFractionCSW ) * ( 9.D0 * gAux1 + &
  &     4.D0 * gAux2 ) + 2.D0 * d2EffPFractionCSW_d2Volume_dDensity * gAux2 * EffPackingFractionCSW * ( 2.D0 + &
  &     EffPackingFractionCSW ) ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) )
  d2AuxEffRDFCSW_d2Volume_dDensity(3) = ( 16.D0 * dxEffPFractionCSW_dVolume_dDensity * dEffPackingFractionCSW_dVolume * &
  &     EffPackingFractionCSW * ( 2.D0 + EffPackingFractionCSW ) * gAux2 + 4.D0 * dEffPackingFractionCSW_dVolume * &
  &     dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dDensity * ( 1.D0 + EffPackingFractionCSW ) * ( 9.D0 * gAux1 + &
  &     8.D0 * gAux2 ) + 8.D0 * dEffPackingFractionCSW_dDensity * d2EffPackingFractionCSW_d2Volume * EffPackingFractionCSW * &
  &     ( 2.D0 + EffPackingFractionCSW ) * gAux2 ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * &
  &     ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) + ( 40.D0 * &
  &     gAux2 * dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dVolume * dEffPackingFractionCSW_dDensity * &
  &     EffPackingFractionCSW * ( 2.D0 + EffPackingFractionCSW ) ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) )
END IF

! First derivative with respect to the volume of the cross derivative of the effective radial distribution function with respect to the volume and the density [mol² . (Å³ / m⁶)]
d2EffectiveRDF_d2Volume_dDensity = SUM( d2AuxEffRDF_d2Volume_dDensity )
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2EffectiveRDFSL_d2Volume_dDensity = SUM( d2AuxEffRDFSutherland_d2Volume_dDensity )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2EffectiveRDFYK_d2Volume_dDensity = SUM( d2AuxEffRDFYukawa_d2Volume_dDensity )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2EffectiveRDFCSW_d2Volume_dDensity = SUM( d2AuxEffRDFCSW_d2Volume_dDensity )
END IF

! First derivative with respect to the volume of the cross derivative of the first-order perturbation coefficient with respect to the volume and the density [K . mol² . (Å³ / m⁶)]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  d2FirstTPTCoeff_d2Volume_dDensity = ( d2cMeanAttractiveEnergy_d2Volume / rNumberDensity ) - ( 2.D0 * &
  &     dcMeanAttractiveEnergy_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) - ( cMeanAttractiveEnergy &
  &     * d2NumberDensity_d2Volume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * cMeanAttractiveEnergy * &
  &     dNumberDensity_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( 1.D0 / &
  &     ( cEffectiveRadialDistributionFunction * cEffectiveRadialDistributionFunction ) ) * dcEffectiveRDF_dVolume * &
  &     ( ( dcMeanAttractiveEnergy_dVolume * dEffRadialDistributionFunct_dDensity ) - ( ( cMeanAttractiveEnergy / &
  &     cEffectiveRadialDistributionFunction ) * dEffRadialDistributionFunct_dDensity * dcEffectiveRDF_dVolume ) + &
  &     ( cMeanAttractiveEnergy * dxEffectiveRDF_dVolume_dDensity ) ) ) + ( ( 1.D0 / cEffectiveRadialDistributionFunction ) * &
  &     ( ( d2cMeanAttractiveEnergy_d2Volume * dEffRadialDistributionFunct_dDensity ) + ( 2.D0 * dcMeanAttractiveEnergy_dVolume * &
  &     dxEffectiveRDF_dVolume_dDensity ) - ( ( 1.D0 / cEffectiveRadialDistributionFunction ) * ( ( dcMeanAttractiveEnergy_dVolume &
  &     * dEffRadialDistributionFunct_dDensity * dcEffectiveRDF_dVolume ) + ( cMeanAttractiveEnergy * dcEffectiveRDF_dVolume * &
  &     dxEffectiveRDF_dVolume_dDensity ) + ( dEffRadialDistributionFunct_dDensity * cMeanAttractiveEnergy * &
  &     d2cEffectiveRDF_d2Volume ) ) ) + ( ( 1.D0 / ( cEffectiveRadialDistributionFunction * cEffectiveRadialDistributionFunction &
  &     ) ) * ( cMeanAttractiveEnergy * dEffRadialDistributionFunct_dDensity * dcEffectiveRDF_dVolume * dcEffectiveRDF_dVolume ) ) &
  &     + ( cMeanAttractiveEnergy * d2EffectiveRDF_d2Volume_dDensity ) ) )
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  d2FirstTPTCoeff_d2Volume_dDensity = ( d2cMeanAttractiveEnergySL_d2Volume / rNumberDensity ) - ( 2.D0 * &
  &     dcMeanAttractiveEnergySutherland_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) - ( &
  &     cMeanAttractiveEnergySutherland * d2NumberDensity_d2Volume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * &
  &     cMeanAttractiveEnergySutherland * dNumberDensity_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity * &
  &     rNumberDensity ) ) - ( ( 1.D0 / ( cEffectiveRadialDistributionFunctionSL * cEffectiveRadialDistributionFunctionSL ) ) * &
  &     dcEffectiveRDFSutherland_dVolume * ( ( dcMeanAttractiveEnergySutherland_dVolume * dEffRadialDistributionFunctSL_dDensity ) &
  &     - ( ( cMeanAttractiveEnergySutherland / cEffectiveRadialDistributionFunctionSL ) * dEffRadialDistributionFunctSL_dDensity &
  &     * dcEffectiveRDFSutherland_dVolume ) + ( cMeanAttractiveEnergySutherland * dxEffectiveRDFSutherland_dVolume_dDensity ) ) ) &
  &     + ( ( 1.D0 / cEffectiveRadialDistributionFunctionSL ) * ( ( d2cMeanAttractiveEnergySL_d2Volume * &
  &     dEffRadialDistributionFunctSL_dDensity ) + ( 2.D0 * dcMeanAttractiveEnergySutherland_dVolume * &
  &     dxEffectiveRDFSutherland_dVolume_dDensity ) - ( ( 1.D0 / cEffectiveRadialDistributionFunctionSL ) * ( ( &
  &     dcMeanAttractiveEnergySutherland_dVolume * dEffRadialDistributionFunctSL_dDensity * dcEffectiveRDFSutherland_dVolume ) + &
  &     ( cMeanAttractiveEnergySutherland * dcEffectiveRDFSutherland_dVolume * dxEffectiveRDFSutherland_dVolume_dDensity ) + &
  &     ( dEffRadialDistributionFunctSL_dDensity * cMeanAttractiveEnergySutherland * d2cEffectiveRDFSutherland_d2Volume ) ) ) + &
  &     ( ( 1.D0 / ( cEffectiveRadialDistributionFunctionSL * cEffectiveRadialDistributionFunctionSL ) ) * &
  &     ( cMeanAttractiveEnergySutherland * dEffRadialDistributionFunctSL_dDensity * dcEffectiveRDFSutherland_dVolume * &
  &     dcEffectiveRDFSutherland_dVolume ) ) + ( cMeanAttractiveEnergySutherland * d2EffectiveRDFSL_d2Volume_dDensity ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  d2FirstTPTCoeff_d2Volume_dDensity = ( d2cMeanAttractiveEnergyYK_d2Volume / rNumberDensity ) - ( 2.D0 * &
  &     dcMeanAttractiveEnergyYukawa_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity ) ) - ( &
  &     cMeanAttractiveEnergyYukawa * d2NumberDensity_d2Volume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * &
  &     cMeanAttractiveEnergyYukawa * dNumberDensity_dVolume * dNumberDensity_dVolume / ( rNumberDensity * rNumberDensity * &
  &     rNumberDensity ) ) - ( ( 1.D0 / ( cEffectiveRadialDistributionFunctionYK * cEffectiveRadialDistributionFunctionYK ) ) * &
  &     dcEffectiveRDFYukawa_dVolume * ( ( dcMeanAttractiveEnergyYukawa_dVolume * dEffRadialDistributionFunctYK_dDensity ) - ( ( &
  &     cMeanAttractiveEnergyYukawa / cEffectiveRadialDistributionFunctionYK ) * dEffRadialDistributionFunctYK_dDensity * &
  &     dcEffectiveRDFYukawa_dVolume ) + ( cMeanAttractiveEnergyYukawa * dxEffectiveRDFYukawa_dVolume_dDensity ) ) ) + ( ( 1.D0 / &
  &     cEffectiveRadialDistributionFunctionYK ) * ( ( d2cMeanAttractiveEnergyYK_d2Volume * dEffRadialDistributionFunctYK_dDensity &
  &     ) + ( 2.D0 * dcMeanAttractiveEnergyYukawa_dVolume * dxEffectiveRDFYukawa_dVolume_dDensity ) - ( ( 1.D0 / &
  &     cEffectiveRadialDistributionFunctionYK ) * ( ( dcMeanAttractiveEnergyYukawa_dVolume * &
  &     dEffRadialDistributionFunctYK_dDensity * dcEffectiveRDFYukawa_dVolume ) + ( cMeanAttractiveEnergyYukawa * &
  &     dcEffectiveRDFYukawa_dVolume * dxEffectiveRDFYukawa_dVolume_dDensity ) + ( dEffRadialDistributionFunctYK_dDensity * &
  &     cMeanAttractiveEnergyYukawa * d2cEffectiveRDFYukawa_d2Volume ) ) ) + ( ( 1.D0 / ( cEffectiveRadialDistributionFunctionYK * &
  &     cEffectiveRadialDistributionFunctionYK ) ) * ( cMeanAttractiveEnergyYukawa * dEffRadialDistributionFunctYK_dDensity * &
  &     dcEffectiveRDFYukawa_dVolume * dcEffectiveRDFYukawa_dVolume ) ) + ( cMeanAttractiveEnergyYukawa * &
  &     d2EffectiveRDFYK_d2Volume_dDensity ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  d2FirstTPTCoeff_d2Volume_dDensity = ( d2cMeanAttractiveEnergyCSW_d2Volume / rNumberDensity ) - ( dNumberDensity_dVolume * &
  &     dcMeanAttractiveEnergyCSW_dVolume / ( rNumberDensity * rNumberDensity ) ) + ( 2.D0 * dNumberDensity_dVolume * &
  &     dNumberDensity_dVolume * cMeanAttractiveEnergyCSW / ( rNumberDensity * rNumberDensity * rNumberDensity ) ) - ( ( &
  &     dcMeanAttractiveEnergyCSW_dVolume * dNumberDensity_dVolume + cMeanAttractiveEnergyCSW * d2NumberDensity_d2Volume ) / ( &
  &     rNumberDensity * rNumberDensity ) ) + ( ( ( d2cMeanAttractiveEnergyCSW_d2Volume * dEffRadialDistributionFunctCSW_dDensity &
  &     ) + ( dcMeanAttractiveEnergyCSW_dVolume * dxEffectiveRDFCSW_dVolume_dDensity ) - ( ( dcMeanAttractiveEnergyCSW_dVolume * &
  &     dcEffectiveRDFCSW_dVolume * dEffRadialDistributionFunctCSW_dDensity / cEffectiveRadialDistributionFunctionCSW ) + &
  &     ( cMeanAttractiveEnergyCSW * d2cEffectiveRDFCSW_d2Volume * dEffRadialDistributionFunctCSW_dDensity / &
  &     cEffectiveRadialDistributionFunctionCSW ) + ( cMeanAttractiveEnergyCSW * dcEffectiveRDFCSW_dVolume * &
  &     dxEffectiveRDFCSW_dVolume_dDensity / cEffectiveRadialDistributionFunctionCSW ) - ( cMeanAttractiveEnergyCSW * &
  &     dcEffectiveRDFCSW_dVolume * dcEffectiveRDFCSW_dVolume * dEffRadialDistributionFunctCSW_dDensity / ( &
  &     cEffectiveRadialDistributionFunctionCSW * cEffectiveRadialDistributionFunctionCSW ) ) ) + &
  &     ( dcMeanAttractiveEnergyCSW_dVolume * dxEffectiveRDFCSW_dVolume_dDensity ) + ( cMeanAttractiveEnergyCSW * &
  &     d2EffectiveRDFCSW_d2Volume_dDensity ) ) / cEffectiveRadialDistributionFunctionCSW ) - ( ( ( ( &
  &     dcMeanAttractiveEnergyCSW_dVolume * dEffRadialDistributionFunctCSW_dDensity ) - ( cMeanAttractiveEnergyCSW * &
  &     dcEffectiveRDFCSW_dVolume * dEffRadialDistributionFunctCSW_dDensity / cEffectiveRadialDistributionFunctionCSW ) + &
  &     ( cMeanAttractiveEnergyCSW * dxEffectiveRDFCSW_dVolume_dDensity ) ) * dcEffectiveRDFCSW_dVolume ) / &
  &     ( cEffectiveRadialDistributionFunctionCSW * cEffectiveRadialDistributionFunctionCSW ) )
END IF

! Second derivative of the mean-attractive energy fluctuation with respect to the volume [K² . mol² / m⁶]
IF( ZhangCorrectionLogical ) THEN
  IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * rNumberDensity * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity * ZhangFactor ) + ( 2.D0 * ZhangFactor * &
    &     dcMeanAttractiveEnergy_dVolume * dHSIsothermalCompressibility_dVolume ) + ( 2.D0 * dHSIsothermalCompressibility_dVolume &
    &     * dZhangFactor_dVolume * dMeanAttractiveEnergy_dDensity * rNumberDensity ) + ( 2.D0 * HSIsothermalCompressibility * &
    &     dZhangFactor_dVolume * dcMeanAttractiveEnergy_dVolume ) + ( 2.D0 * dxFirstTPTCoefficient_dVolume_dDensity * &
    &     rNumberDensity * dZhangFactor_dVolume * HSIsothermalCompressibility ) + ( HSIsothermalCompressibility * &
    &     dxFirstTPTCoefficient_dVolume_dDensity * ZhangFactor * dNumberDensity_dVolume ) + ( HSIsothermalCompressibility * &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * d2ZhangFactor_d2Volume ) + ( d2HSIsothermalCompressibility_d2Volume * &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * ZhangFactor ) + ( d2cMeanAttractiveEnergy_d2Volume * &
    &     HSIsothermalCompressibility * ZhangFactor ) + ( rNumberDensity * HSIsothermalCompressibility * &
    &     d2FirstTPTCoeff_d2Volume_dDensity * ZhangFactor ) )
  ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity * ZhangFactor ) + ( 2.D0 &
    &     * dHSIsothermalCompressibility_dVolume * dcMeanAttractiveEnergySutherland_dVolume * ZhangFactor ) + ( 2.D0 * &
    &     rNumberDensity * dHSIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity * dZhangFactor_dVolume ) + ( 2.D0 &
    &     * dZhangFactor_dVolume * HSIsothermalCompressibility * dcMeanAttractiveEnergySutherland_dVolume ) + ( 2.D0 * &
    &     HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity * dZhangFactor_dVolume ) + ( &
    &     dxFirstTPTCoefficient_dVolume_dDensity * HSIsothermalCompressibility * ZhangFactor * dNumberDensity_dVolume ) + ( &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * d2ZhangFactor_d2Volume * HSIsothermalCompressibility ) + ( &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * d2HSIsothermalCompressibility_d2Volume * ZhangFactor ) + ( &
    &     HSIsothermalCompressibility * d2cMeanAttractiveEnergySL_d2Volume * ZhangFactor ) + ( rNumberDensity * &
    &     HSIsothermalCompressibility * d2FirstTPTCoeff_d2Volume_dDensity * ZhangFactor ) )
  ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity * ZhangFactor ) + ( 2.D0 &
    &     * dHSIsothermalCompressibility_dVolume * dcMeanAttractiveEnergyYukawa_dVolume * ZhangFactor ) + ( 2.D0 * rNumberDensity &
    &     * dHSIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity * dZhangFactor_dVolume ) + ( 2.D0 * &
    &     dZhangFactor_dVolume * HSIsothermalCompressibility * dcMeanAttractiveEnergyYukawa_dVolume ) + ( 2.D0 * &
    &     HSIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity * dZhangFactor_dVolume ) + ( &
    &     dxFirstTPTCoefficient_dVolume_dDensity * HSIsothermalCompressibility * ZhangFactor * dNumberDensity_dVolume ) + ( &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * d2ZhangFactor_d2Volume * HSIsothermalCompressibility ) + ( &
    &     rNumberDensity * dMeanAttractiveEnergy_dDensity * d2HSIsothermalCompressibility_d2Volume * ZhangFactor ) + ( &
    &     HSIsothermalCompressibility * d2cMeanAttractiveEnergyYK_d2Volume * ZhangFactor ) + ( rNumberDensity * &
    &     HSIsothermalCompressibility * d2FirstTPTCoeff_d2Volume_dDensity * ZhangFactor ) )
  ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * rNumberDensity * &
    &     dxFirstTPTCoefficient_dVolume_dDensity * dHCBIsothermalCompressibility_dVolume * ZhangFactor ) + ( 2.D0 * &
    &     dcMeanAttractiveEnergyCSW_dVolume * dHCBIsothermalCompressibility_dVolume * ZhangFactor ) + &
    &     ( dMeanAttractiveEnergy_dDensity * d2HCBIsothermalCompressibility_d2Volume * rNumberDensity * ZhangFactor ) + &
    &     ( d2cMeanAttractiveEnergyCSW_d2Volume * HCBIsothermalCompressibility * ZhangFactor ) + ( HCBIsothermalCompressibility * &
    &     rNumberDensity * d2FirstTPTCoeff_d2Volume_dDensity * ZhangFactor ) + ( HCBIsothermalCompressibility * &
    &     dxFirstTPTCoefficient_dVolume_dDensity * dNumberDensity_dVolume * ZhangFactor ) + ( 2.D0 * rNumberDensity * &
    &     dHCBIsothermalCompressibility_dVolume * dMeanAttractiveEnergy_dDensity * dZhangFactor_dVolume ) + ( 2.D0 * &
    &     dZhangFactor_dVolume * HCBIsothermalCompressibility * dcMeanAttractiveEnergyCSW_dVolume ) + ( 2.D0 * &
    &     HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity * dZhangFactor_dVolume ) + &
    &     ( rNumberDensity * dMeanAttractiveEnergy_dDensity * d2ZhangFactor_d2Volume * HCBIsothermalCompressibility ) )
  END IF
ELSE
  IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * rNumberDensity * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity ) + ( 2.D0 * &
    &     dcMeanAttractiveEnergy_dVolume * dHSIsothermalCompressibility_dVolume ) + ( HSIsothermalCompressibility * &
    &     dxFirstTPTCoefficient_dVolume_dDensity * dNumberDensity_dVolume ) + ( rNumberDensity * &
    &     d2HSIsothermalCompressibility_d2Volume * dMeanAttractiveEnergy_dDensity ) + ( d2cMeanAttractiveEnergy_d2Volume * &
    &     HSIsothermalCompressibility ) + ( rNumberDensity * HSIsothermalCompressibility * d2FirstTPTCoeff_d2Volume_dDensity ) )
  ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity ) + ( 2.D0 * &
    &     dcMeanAttractiveEnergySutherland_dVolume * dHSIsothermalCompressibility_dVolume ) + ( &
    &     dxFirstTPTCoefficient_dVolume_dDensity * HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * &
    &     d2HSIsothermalCompressibility_d2Volume * dMeanAttractiveEnergy_dDensity ) + ( d2cMeanAttractiveEnergySL_d2Volume * &
    &     HSIsothermalCompressibility ) + ( rNumberDensity * d2FirstTPTCoeff_d2Volume_dDensity * HSIsothermalCompressibility ) )
  ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * &
    &     dHSIsothermalCompressibility_dVolume * dxFirstTPTCoefficient_dVolume_dDensity * rNumberDensity ) + ( 2.D0 * &
    &     dcMeanAttractiveEnergyYukawa_dVolume * dHSIsothermalCompressibility_dVolume ) + ( dxFirstTPTCoefficient_dVolume_dDensity &
    &     * HSIsothermalCompressibility * dNumberDensity_dVolume ) + ( rNumberDensity * d2HSIsothermalCompressibility_d2Volume * &
    &     dMeanAttractiveEnergy_dDensity ) + ( d2cMeanAttractiveEnergyYK_d2Volume * HSIsothermalCompressibility ) + ( &
    &     rNumberDensity * d2FirstTPTCoeff_d2Volume_dDensity * HSIsothermalCompressibility ) )
  ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
    d2cMeanEnergyFluctuations_d2Volume = 0.5D0 * ijaWellDepth(cComponent,cComponent) * ( ( 2.D0 * rNumberDensity * &
    &     dxFirstTPTCoefficient_dVolume_dDensity * dHCBIsothermalCompressibility_dVolume ) + ( 2.D0 * &
    &     dcMeanAttractiveEnergyCSW_dVolume * dHCBIsothermalCompressibility_dVolume ) + ( dMeanAttractiveEnergy_dDensity * &
    &     d2HCBIsothermalCompressibility_d2Volume * rNumberDensity ) + ( d2cMeanAttractiveEnergyCSW_d2Volume * &
    &     HCBIsothermalCompressibility ) + ( HCBIsothermalCompressibility * rNumberDensity * d2FirstTPTCoeff_d2Volume_dDensity ) + &
    &     ( HCBIsothermalCompressibility * dxFirstTPTCoefficient_dVolume_dDensity * dNumberDensity_dVolume ) )
  END IF
END IF

! Second derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MeanFluctuationFreeEnergy_d2Volume = d2cMeanEnergyFluctuations_d2Volume * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
IF( HigherOrderTPTLogical ) THEN
  d2HigherOrderTerms_d2Volume = 0.D0
  DO nOrder = 3, nHigherOrder
    d2HigherOrderTerms_d2Volume = d2HigherOrderTerms_d2Volume + ( d2MeanAttractiveFreeEnergy_d2Volume * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) ) + ( dMeanAttractiveFreeEnergy_dVolume * ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * &
    &     ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy &
    &     * dMeanAttractiveFreeEnergy_dVolume) ) / ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) &
    &     + ( ( 4.D0 * DBLE( nOrder - 1 ) * DBLE( nOrder - 2 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 3 ) ) * &
    &     ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - &
    &     (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) * ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) ) / &
    &     ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) + &
    &     ( ( 2.D0 * DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / &
    &     mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     d2MeanFluctuationFreeEnergy_d2Volume) - (mMeanAttFluctuationHelmholtzFreeEnergy * d2MeanAttractiveFreeEnergy_d2Volume) ) &
    &     / ( mMeanAttractiveHelmholtzFreeEnergy ) ) - ( ( ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dMeanFluctuationFreeEnergy_dVolume) - (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) * &
    &     dMeanAttractiveFreeEnergy_dVolume ) / ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) ) )
  END DO
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
  d2HigherOrderFEnergy_d2Volume = d2HigherOrderTerms_d2Volume ! Proven units
ELSE
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
  d2HigherOrderFEnergy_d2Volume = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Second derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2MonomerFreeEnergy_d2Volume = d2ExcludedVolumeFreeEnergy_d2Volume + d2MeanAttractiveFreeEnergy_d2Volume + &
&     d2MeanFluctuationFreeEnergy_d2Volume + d2HigherOrderFEnergy_d2Volume ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Second derivative of the total Helmholtz free energy with respect to the volume [Pa . mol / m³]
d2TotalFreeEnergy_d2Volume = d2IdealFreeEnergy_d2Volume + d2MonomerFreeEnergy_d2Volume ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!  CROSS DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE AND THE VOLUME   !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Cross derivative of the ideal contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxIdealFreeEnergy_dVolume_dTemperature = dIdealFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Cross derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxEVFreeEnergy_dVolume_dTemperature = dExcludedVolumeFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Cross derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxSecondOrderTPT_dVolume_dTemperature = - dMeanFluctuationFreeEnergy_dVolume / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
IF( HigherOrderTPTLogical ) THEN
  dxHigherOrderTerms_dVolume_dTemperature = 0.D0
  DO nOrder = 3, nHigherOrder
    ! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
    dMeanFluctuationFreeEnergy_dTemperature = - cMeanAttractiveEnergyFluctuations * cUniversalGas / ( Temperature * Temperature ) ! Proven units
    ! Barker-Henderson's higher-order perturbation theory (approximation)
    dxHigherOrderTerms_dVolume_dTemperature = dxHigherOrderTerms_dVolume_dTemperature + ( 2.D0 * dMeanAttractiveFreeEnergy_dVolume &
    &     * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / &
    &     mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * (dMeanFluctuationFreeEnergy_dTemperature / &
    &     mMeanAttractiveHelmholtzFreeEnergy) ) + ( 4.D0 * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 3 ) ) / &
    &     ( mMeanAttractiveHelmholtzFreeEnergy * mMeanAttractiveHelmholtzFreeEnergy ) ) * ( ( DBLE( nOrder - 2 ) * &
    &     dMeanFluctuationFreeEnergy_dTemperature * ( (mMeanAttractiveHelmholtzFreeEnergy * dMeanFluctuationFreeEnergy_dVolume) - &
    &     (mMeanAttFluctuationHelmholtzFreeEnergy * dMeanAttractiveFreeEnergy_dVolume) ) ) + &
    &     ( mMeanAttFluctuationHelmholtzFreeEnergy * ( (mMeanAttractiveHelmholtzFreeEnergy * &
    &     dxSecondOrderTPT_dVolume_dTemperature) - (dMeanFluctuationFreeEnergy_dTemperature * &
    &     dMeanAttractiveFreeEnergy_dVolume) ) ) )
  END DO
  ! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
  dxHigherOrderFEnergy_dVolume_dTemperature = dxHigherOrderTerms_dVolume_dTemperature ! Proven units
ELSE
  ! Cross derivative of the Barker-Henderson's higher-order perturbation terms with respect to the volume and the temperature [Pa / K]
  dxHigherOrderFEnergy_dVolume_dTemperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Cross derivative of the total monomer-monomer contribution to the Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxMonomerFEnergy_dVolume_dTemperature = dxEVFreeEnergy_dVolume_dTemperature + dxSecondOrderTPT_dVolume_dTemperature + &
&     dxHigherOrderFEnergy_dVolume_dTemperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Cross derivative of the total Helmholtz free energy with respect to the volume and the temperature [Pa / K]
dxTotalFreeEnergy_dVolume_dTemperature = dxIdealFreeEnergy_dVolume_dTemperature + dxMonomerFEnergy_dVolume_dTemperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                     THERMODYNAMIC PROPERTIES                                     !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Isothermal compressibility [1 / Pa]
IsothermalCompressibility = 1.D0 / ( mVolume * d2TotalFreeEnergy_d2Volume )

! Thermal expansion coefficient [1 / K]
ThermalExpansionCoefficient = - IsothermalCompressibility * dxTotalFreeEnergy_dVolume_dTemperature

! Pressure [Pa]
Pressure = - dTotalFreeEnergy_dVolume

! Compressibility factor [unitless]
CompressibilityFactor = Pressure * mVolume / ( cUniversalGas * Temperature )

RETURN

END SUBROUTINE Calculate_Pressure_Single_Component