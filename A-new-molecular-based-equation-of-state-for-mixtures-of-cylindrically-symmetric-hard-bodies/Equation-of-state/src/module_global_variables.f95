! ************************************************************************************************ !
!                                         GLOBAL VARIABLES                                         !
! ************************************************************************************************ !
!       This module defines a few global variables, strings, and constants to be used in the       !
!                           main program and its underlying subroutines.                           !
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
MODULE GlobalVariables

! Uses an intrinsic module to define the kind of variables
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Int64, Real64

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: nComponents            ! Number of components
INTEGER( Kind= Int64 ) :: PotentialType          ! Potential type (1: Square-Well, 2: Sutherland, 3: Yukawa)
INTEGER( Kind= Int64 ) :: nSimplexCycles         ! Number of simplex cycles
INTEGER( Kind= Int64 ) :: NonSphericalMixingRule ! Mixing rule for non-spherical elements (1: mixing rule applied directly to the aspect ratio, 2: mixing rule applied directly to the angle average of the excluded volume of a pair of nonspherical particles)
INTEGER( Kind= Int64 ) :: EffPFractionMixingRule ! Mixing rule for the effective packing fraction (1: uses the reduced density 3 of the mixture of spherical segments, 2: uses the one-fluid van der Waals mixing rule)
INTEGER( Kind= Int64 ) :: nHigherOrder           ! Number of high-order perturbation terms (starting from the third TPT coefficient)

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: ijWellDepthCorrection      ! Correction to the well depth of the interaction between components
REAL( Kind= Real64 ) :: ijPotentialRangeCorrection ! Correction to the potential range of the interaction between components
REAL( Kind= Real64 ) :: LastObjectiveFunction      ! Last value of the objective function

! ************************************************************************************************ !
! REAL VARIABLES (ALLOCATABLE)                                                                     !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: TemperatureParameter ! Temperature-dependent coefficients (NIST or TRC)

! ************************************************************************************************ !
! REAL VARIABLES (ALLOCATABLE)                                                                     !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cGeometry                      ! Molecular geometry of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cMolarMass                     ! Molar mass of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cAspectRatio                   ! Aspect ratio of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: fAspectRatio                   ! Aspect ratio of a component field
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cDiameter                      ! Molecular diameter of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: fDiameter                      ! Molecular diameter of a component field
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cDiameterSphere                ! Molecular diameter of a sphere with the same volume of the non-spherical geometry
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aDiameter                      ! Molecular diameter of a component [Å]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aDiameterField                 ! Molecular diameter of a component field [Å]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aDiameterSphere                ! Molecular diameter of a sphere with the same volume of the non-spherical geometry [Å]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cLength                        ! Molecular length of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: fLength                        ! Molecular length of a component field
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aLength                        ! Molecular length of a component [Å]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aLengthField                   ! Molecular length of a component field [Å]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cWellDepth                     ! Well depth of a component
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: aWellDepth                     ! Well depth of a component [K]
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: cPotentialRange                ! Potential range of a component
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijAspectRatio                  ! Aspect ratio of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijDiameter                     ! Diameter of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijaDiameter                    ! Diameter of the interaction between components (reduced)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijDiameterSphere               ! Spherical diameter of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijaDiameterSphere              ! Spherical diameter of the interaction between components (reduced)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijDiameterSphereCubic          ! Spherical dimater of the interaction between components (cubic)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijaDiameterSphereCubic         ! Spherical diameter of the interaction between components (cubic, reduced)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijLength                       ! Length of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijaLength                      ! Length of the interaction between components (reduced)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijPotentialRange               ! Range of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijPotentialRangeSquared        ! Range of the interaction between components (squared)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijPotentialRangeCubic          ! Range of the interaction between components (cubic)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijWellDepth                    ! Well depth of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijaWellDepth                   ! Well depth of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijSecondVirialCoefficient      ! Second virial coefficient of the interaction between components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijSecondVirialCoefficientField ! Second virial coefficient of the interaction between components (convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijNonSphericity                ! Nonsphericity parameter of the interaction between components (convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijHSSecondVirialCoefficient    ! Second virial coefficient of the interaction between equivalent hard spherical components
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ijRatioSecondVirialCoefficient ! Ratio of second virial coefficients of the interaction between components

! ************************************************************************************************ !
! REAL PARAMETERS (CONSTANTS)                                                                      !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: cBoltzmann    = 1.38064852D-23         ! Boltzmann's constant in Joule / Kelvin
REAL( Kind= Real64 ), PARAMETER :: cAvogadro     = 6.02214076D+23         ! Avogadro's constant in 1 / mol
REAL( Kind= Real64 ), PARAMETER :: cPlanck       = 6.62607004D-34         ! Planck's constant in Joule . second
REAL( Kind= Real64 ), PARAMETER :: cPi           = 4.D0 * DATAN( 1.D0 )   ! π number
REAL( Kind= Real64 ), PARAMETER :: cUniversalGas = cAvogadro * cBoltzmann ! Universal gas constant in Joule / (Kelvin . mol)

! ************************************************************************************************ !
! REAL PARAMETERS (OTHER CONSTANTS)                                                                !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: cDeBroglie      = cPlanck * cAvogadro         ! Constant in (Joule / mol) . second
REAL( Kind= Real64 ), PARAMETER :: cDeBroglieCubic = cDeBroglie ** 3.D0          ! Constant in [(Joule / mol) . second]³
REAL( Kind= Real64 ), PARAMETER :: cDeBroglieIdeal = cDeBroglieCubic * cAvogadro ! Constant in [(Joule³ / mol⁴) . second³]

! ************************************************************************************************ !
! REAL PARAMETERS (CSW)                                                                            !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWEORCoefficientA1 = RESHAPE( [ [2.258550D0, 10.843815D0, -5.797034D0], &
&                                                                                      [-1.50349D0, -13.76179D0, 7.2894411D0], &
&                                                                                      [0.249434D0, 4.3221109D0, -2.270221D0], &
&                                                                                      [-0.669270D0, -65.05305D0, 25.258131D0], &
&                                                                                      [1.4004900D0, 81.934716D0, -31.61837D0], &
&                                                                                      [-0.827739D0, -25.44501D0, 9.7403240D0], &
&                                                                                      [1.015760D1, 82.813878D0, -5.193670D0], &
&                                                                                      [-1.50427D1, -105.3865D0, 5.1256200D0], &
&                                                                                      [5.308270D0, 33.020905D0, -1.295402D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction (ellipsoids-of-revolution)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWEORCoefficientA2 = RESHAPE( [ [8.3786934D0, 11.834703D0, -4.919689D0], &
&                                                                                      [-7.477416D0, -14.71838D0, 6.1138737D0], &
&                                                                                      [1.5175421D0, 4.5985847D0, -1.835521D0], &
&                                                                                      [-4.560861D0, -72.34344D0, 22.352753D0], &
&                                                                                      [1.7964686D0, 90.818316D0, -27.75594D0], &
&                                                                                      [0.5040354D0, -27.31398D0, 7.9427089D0], &
&                                                                                      [-0.874200D0, 80.621100D0, -7.211933D0], &
&                                                                                      [1.9707776D0, -100.7248D0, 8.2342050D0], &
&                                                                                      [-0.282994D0, 29.796480D0, -1.323709D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction without Percus-Yevick correction to the isothermal compressibility of HCBs (ellipsoids-of-revolution)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWEORCoefficientA2PY = RESHAPE( [ [8.3401990D0, 11.702345D0, -5.028671D0], &
&                                                                                        [-7.418345D0, -14.57155D0, 6.2736036D0], &
&                                                                                        [1.4919418D0, 4.5528857D0, -1.882892D0], &
&                                                                                        [-4.168421D0, -72.01220D0, 22.909173D0], &
&                                                                                        [1.3271749D0, 90.607899D0, -28.67295D0], &
&                                                                                        [0.6776622D0, -27.23982D0, 8.2248570D0], &
&                                                                                        [-1.633309D0, 80.115614D0, -6.928333D0], &
&                                                                                        [2.9172040D0, -100.3804D0, 8.3392552D0], &
&                                                                                        [-0.578872D0, 29.665245D0, -1.383627D0]  &
&                                                                                      ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction with Percus-Yevick correction to the isothermal compressibility of HCBs (ellipsoids-of-revolution)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWSPCCoefficientA1 = RESHAPE( [ [1.6929073D0, 11.164094D0, -5.8771841D0], &
&                                                                                      [-0.793720D0, -14.60521D0, 7.62976630D0], &
&                                                                                      [0.0337918D0, 4.6767380D0, -2.4303734D0], &
&                                                                                      [5.2858575D0, -62.20125D0, 23.3968270D0], &
&                                                                                      [-6.244525D0, 80.072894D0, -30.161203D0], &
&                                                                                      [1.5720075D0, -25.41248D0, 9.57507900D0], &
&                                                                                      [2.0705564D0, 79.325009D0, -7.6908206D0], &
&                                                                                      [-4.984929D0, -103.7028D0, 9.52653130D0], &
&                                                                                      [2.2633304D0, 33.443747D0, -3.0599829D0]  &
&                                                                                    ], [3, 9] )                              ! Matrix of coefficients of the effective packing fraction (spherocylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWSPCCoefficientA2 = RESHAPE( [ [8.5341889D0, 12.799608D0, -4.9822504D0], &
&                                                                                      [-7.697638D0, -16.04470D0, 6.23838980D0], &
&                                                                                      [1.5831103D0, 5.0398762D0, -1.9006712D0], &
&                                                                                      [-4.735795D0, -73.30082D0, 20.8021480D0], &
&                                                                                      [2.1439287D0, 91.013795D0, -25.715493D0], &
&                                                                                      [0.4175330D0, -27.36730D0, 7.48010160D0], &
&                                                                                      [-1.409160D0, 80.970946D0, -7.5175184D0], &
&                                                                                      [2.5296130D0, -99.51427D0, 9.21575340D0], &
&                                                                                      [-0.495048D0, 29.440899D0, -2.2100362D0]  &
&                                                                                    ], [3, 9] )                              ! Matrix of coefficients of the effective packing fraction without Percus-Yevick correction to the isothermal compressibility of HCBs (spherocylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWSPCCoefficientA2PY = RESHAPE( [ [ 8.287197D0, 12.923605D0, -5.032596D0], &
&                                                                                        [-7.366823D0, -16.32380D0, 6.3641360D0], &
&                                                                                        [ 1.474453D0, 5.1531012D0, -1.953500D0], &
&                                                                                        [-3.632363D0, -72.89690D0, 21.090014D0], &
&                                                                                        [ 0.725000D0, 91.180185D0, -26.32737D0], &
&                                                                                        [ 0.885557D0, -27.54556D0, 7.7165702D0], &
&                                                                                        [-2.585457D0, 79.355026D0, -8.376672D0], &
&                                                                                        [ 4.080708D0, -98.18759D0, 10.412813D0], &
&                                                                                        [-0.985689D0, 29.140912D0, -2.581224D0]  &
&                                                                                      ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction with Percus-Yevick correction to the isothermal compressibility of HCBs (spherocylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWCYLCoefficientA1 = RESHAPE( [ [1.7364071D0, 7.302352D0, -2.7952615D0], &
&                                                                                      [-0.653163D0, -10.7583D0, 4.43013920D0], &
&                                                                                      [-0.091113D0, 3.730956D0, -1.6018576D0], &
&                                                                                      [6.2750733D0, -63.3687D0, 21.1494450D0], &
&                                                                                      [-9.686810D0, 81.22186D0, -27.499447D0], &
&                                                                                      [3.2605216D0, -25.7433D0, 8.89311610D0], &
&                                                                                      [-3.774129D0, 83.64615D0, -12.432418D0], &
&                                                                                      [3.0687882D0, -98.6021D0, 11.3854950D0], &
&                                                                                      [-0.507658D0, 29.61948D0, -3.0145526D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction (cylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWCYLCoefficientA2 = RESHAPE( [ [7.6875099D0, 7.2650887D0, -1.073949D0], &
&                                                                                      [-6.364866D0, -14.39882D0, 4.8186344D0], &
&                                                                                      [1.0600282D0, 5.5289221D0, -2.108624D0], &
&                                                                                      [-3.945195D0, -69.91426D0, 17.410788D0], &
&                                                                                      [-0.464489D0, 97.243499D0, -28.02057D0], &
&                                                                                      [1.6645452D0, -30.60729D0, 9.0702075D0], &
&                                                                                      [-3.527090D0, 80.905112D0, -9.395649D0], &
&                                                                                      [5.9845264D0, -95.80396D0, 8.9825814D0], &
&                                                                                      [-1.576637D0, 25.950589D0, -0.784142D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction without Percus-Yevick correction to the isothermal compressibility of HCBs (cylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWCYLCoefficientA2PY = RESHAPE( [ [7.6489305D0, 6.9565195D0, -0.908378D0], &
&                                                                                        [-6.194723D0, -14.42684D0, 4.8798028D0], &
&                                                                                        [0.9777838D0, 5.6474343D0, -2.197378D0], &
&                                                                                        [-4.021795D0, -67.45465D0, 16.329071D0], &
&                                                                                        [-1.082980D0, 96.249628D0, -28.12525D0], &
&                                                                                        [2.0347838D0, -30.94730D0, 9.5271617D0], &
&                                                                                        [-3.464830D0, 78.167624D0, -8.788455D0], &
&                                                                                        [6.8706089D0, -95.52149D0, 10.513912D0], &
&                                                                                        [-2.144623D0, 26.866991D0, -1.933863D0]  &
&                                                                                      ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction with Percus-Yevick correction to the isothermal compressibility of HCBs (cylinders)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWAllCoefficientA1 = RESHAPE( [ [1.2861436D0, 8.6662363D0, -3.640347D0], &
&                                                                                      [-0.264516D0, -10.85520D0, 4.3270990D0], &
&                                                                                      [-0.142175D0, 3.4026527D0, -1.300220D0], &
&                                                                                      [8.0689592D0, -65.54128D0, 23.656195D0], &
&                                                                                      [-9.789071D0, 78.593035D0, -26.74209D0], &
&                                                                                      [2.7112206D0, -23.79707D0, 7.6376301D0], &
&                                                                                      [-3.591238D0, 81.118755D0, -12.30145D0], &
&                                                                                      [2.5508089D0, -99.59120D0, 11.502932D0], &
&                                                                                      [-0.232461D0, 30.926322D0, -2.627241D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction for all geometries (Convex square-well potential) | <This includes all data from EOR, SPC and CYL (disks and rods)>
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWAllCoefficientA2 = RESHAPE( [ [7.8994465D0, 13.834713D0, -5.695347D0], &
&                                                                                      [-6.917255D0, -18.70872D0, 7.7527537D0], &
&                                                                                      [1.3686527D0, 6.3188469D0, -2.559332D0], &
&                                                                                      [-1.282246D0, -68.59763D0, 19.349581D0], &
&                                                                                      [-2.259924D0, 92.038508D0, -25.60741D0], &
&                                                                                      [1.7079665D0, -29.99074D0, 7.9999512D0], &
&                                                                                      [-4.140077D0, 75.742761D0, -7.464172D0], &
&                                                                                      [5.8344133D0, -102.0788D0, 9.1983097D0], &
&                                                                                      [-1.378528D0, 33.420542D0, -2.367931D0]  &
&                                                                                    ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction for all geometries without Percus-Yevick correction to the isothermal compressibility of HCBs (Convex square-well potential) | <This includes all data from EOR and SPC>
REAL( Kind= Real64 ), DIMENSION( 3, 9 ), PARAMETER :: CSWAllCoefficientA2PY = RESHAPE( [ [7.6365439D0, 14.255371D0, -6.460901D0], &
&                                                                                        [-6.611867D0, -18.94982D0, 8.5707438D0], &
&                                                                                        [1.2811614D0, 6.2830786D0, -2.747545D0], &
&                                                                                        [0.0310147D0, -67.69533D0, 22.197234D0], &
&                                                                                        [-3.651557D0, 88.699091D0, -27.97630D0], &
&                                                                                        [2.0782421D0, -28.12528D0, 8.2298568D0], &
&                                                                                        [-5.728935D0, 71.519280D0, -10.13859D0], &
&                                                                                        [7.4793834D0, -93.25401D0, 10.594890D0], &
&                                                                                        [-1.763158D0, 29.383458D0, -1.999864D0]  &
&                                                                                      ], [3, 9] )                             ! Matrix of coefficients of the effective packing fraction for all geometries with Percus-Yevick correction to the isothermal compressibility of HCBs (Convex square-well potential) | <This includes all data from EOR and SPC>

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 04 ) :: SpecificHeatReference ! Specific heat reference (TRC or NIST)
CHARACTER( Len= 08 ) :: DescriptorDate        ! Descriptor date
CHARACTER( Len= 08 ) :: DescriptorHour        ! Descriptor hour
CHARACTER( Len= 12 ) :: FormatDate            ! Format date
CHARACTER( Len= 12 ) :: FormatHour            ! Format hour

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 25 ), DIMENSION( : ), ALLOCATABLE :: cMoleculeName ! Name (or formula) of a component
CHARACTER( Len= 10 ), DIMENSION( : ), ALLOCATABLE :: cFormulaName  ! Structural formula of a component
CHARACTER( Len= 23 ), DIMENSION( : ), ALLOCATABLE :: cGeometryName ! Geometry name of a component

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: HigherOrderTPTLogical   ! Use or neglect higher-order terms in the TPT expansion
LOGICAL                 :: ZhangCorrectionLogical  ! Use or neglect Zhang's correction on the second-order perturbation coefficient
LOGICAL                 :: PYHCBCorrectionLogical  ! Use or neglect Percus-Yevick correction to the isothermal compressibility of HCBs
LOGICAL                 :: ReferenceBoublikLogical ! Use or neglect the Boublik expression for the reference system of non-spherical particles
LOGICAL                 :: UseA1ForA2Logical       ! Replace [TRUE] the effective packing fraction coefficients for the second-order perturbation coefficient (CSW potential) by those for the first-order perturbation coefficient (CSW potential)
LOGICAL                 :: UseA1AllGeometries      ! Replace [TRUE] the effective packing fraction coefficients for the first-order perturbation coefficient (CSW potential) by those for all geometries (CSW potential)
LOGICAL                 :: UseA2AllGeometries      ! Replace [TRUE] the effective packing fraction coefficients for the second-order perturbation coefficient (CSW potential) by those for all geometries (CSW potential)
LOGICAL, DIMENSION( 4 ) :: PotentialTypeLogical    ! Logical array for the potential type (1: Spherical Square-Well, 2: Sutherland, 3: Yukawa, 4: Convex Square-Well)

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL, DIMENSION( :, : ), ALLOCATABLE :: GeometrySpecification ! Geometry specification (ellipsoids-of-revolution, spherocylinders or cylinders)

END MODULE GlobalVariables