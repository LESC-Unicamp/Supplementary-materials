! ************************************************************************************************ !
!                                         GENERAL ROUTINES                                         !
! ************************************************************************************************ !
!             This subroutine contains all general routines used in the main program.              !
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
!                  This function converts any string into uppercase (from A to Z)                  !
! ************************************************************************************************ !
SUBROUTINE ToUpper( StringInput, StringLength, StringOutput )

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER :: StringLength ! String length
INTEGER :: iCharacter   ! ASCII character code
INTEGER :: iString      ! Counter

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( LEN= StringLength ) :: StringInput  ! Length of input string
CHARACTER( LEN= StringLength ) :: StringOutput ! Length of output string

! Character positions
DO iString = 1, StringLength
  ! ASCII character code
  iCharacter = IACHAR( StringInput(iString:iString) )
  ! Convert to uppercase (letters only)
  IF( iCharacter >= IACHAR( "a" ) .AND. iCharacter <= IACHAR( "z" ) ) THEN
    StringOutput(iString:iString) = ACHAR(IACHAR(StringInput(iString:iString))-32)
  ! Do not convert symbols or numbers (special characters included)
  ELSE
    StringOutput(iString:iString) = StringInput(iString:iString)
  END IF
END DO

RETURN

END SUBROUTINE ToUpper

! ************************************************************************************************ !
!                  This subroutine allocates the arrays used in the main program                   !
! ************************************************************************************************ !
SUBROUTINE Allocation_Arrays(  )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! Allocation
ALLOCATE( cMoleculeName(nComponents), cFormulaName(nComponents), cMolarMass(nComponents) )
ALLOCATE( cDiameter(nComponents), fDiameter(nComponents), aDiameter(nComponents), aDiameterField(nComponents) )
ALLOCATE( cLength(nComponents), fLength(nComponents), aLength(nComponents), aLengthField(nComponents) )
ALLOCATE( ijDiameter(nComponents,nComponents), ijaDiameter(nComponents,nComponents) )
ALLOCATE( ijLength(nComponents,nComponents), ijaLength(nComponents,nComponents) )
ALLOCATE( cAspectRatio(nComponents), fAspectRatio(nComponents), ijAspectRatio(nComponents,nComponents) )
ALLOCATE( cWellDepth(nComponents), aWellDepth(nComponents) )
ALLOCATE( ijWellDepth(nComponents,nComponents), ijaWellDepth(nComponents,nComponents) )
ALLOCATE( cDiameterSphere(nComponents), aDiameterSphere(nComponents) )
ALLOCATE( cPotentialRange(nComponents), ijPotentialRange(nComponents,nComponents) )
ALLOCATE( ijDiameterSphere(nComponents,nComponents), ijaDiameterSphere(nComponents,nComponents) )
ALLOCATE( ijDiameterSphereCubic(nComponents,nComponents) )
ALLOCATE( ijPotentialRangeSquared(nComponents,nComponents), ijPotentialRangeCubic(nComponents,nComponents) )
ALLOCATE( TemperatureParameter(nComponents,5) )
ALLOCATE( ijSecondVirialCoefficient(nComponents,nComponents), ijHSSecondVirialCoefficient(nComponents,nComponents) )
ALLOCATE( ijSecondVirialCoefficientField(nComponents,nComponents), ijNonSphericity(nComponents,nComponents) )
ALLOCATE( ijRatioSecondVirialCoefficient(nComponents,nComponents) )
ALLOCATE( GeometrySpecification(nComponents,3) )

RETURN

END SUBROUTINE Allocation_Arrays

! ************************************************************************************************ !
!                 This subroutine reads the input file and stores the information                  !
! ************************************************************************************************ !
SUBROUTINE Read_Global_Specification(  )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 01 ) :: Dummy ! Dummy

! Read specifications
OPEN( Unit= 10, File= "specs.ini", Action= "READ" )
READ( 10, * ) Dummy, nComponents
IF( nComponents < 1 ) THEN
  WRITE( *, "(G0)" ) "Error: Number of components must be greater than zero. Exiting..."
  CALL EXIT(  )
END IF
ALLOCATE( cGeometry(nComponents) )
ALLOCATE( cGeometryName(nComponents) )
READ( 10, * ) Dummy, cGeometry(1:nComponents)
IF( ANY( cGeometry < 1 ) .OR. ANY( cGeometry > 3 ) ) THEN
  WRITE( *, "(G0)" ) "Error: Geometry must be 1 (Elliposids-of-Revolution), 2 (Spherocylinders) or 3 (Cylinders). Exiting..."
  CALL EXIT(  )
END IF
DO cComponent = 1, nComponents
  IF( cGeometry(cComponent) == 1 ) THEN
    cGeometryName(cComponent) = "Ellipsoid-of-Revolution"
  ELSE IF( cGeometry(cComponent) == 2 ) THEN
    cGeometryName(cComponent) = "Spherocylinder"
  ELSE IF( cGeometry(cComponent) == 3 ) THEN
    cGeometryName(cComponent) = "Cylinder"
  END IF
END DO
READ( 10, * ) Dummy, NonSphericalMixingRule
IF( NonSphericalMixingRule < 1 .OR. NonSphericalMixingRule > 4 ) THEN
  WRITE( *, "(2G0)" ) "Error: Non-spherical mixing rule must be 1 (aspect ratio rule), 2 (angle average rule), 3 ", &
  &                   "(ratio of second virial coefficients), or 4 (second virial coefficients). Exiting..."
  CALL EXIT(  )
END IF
DO cComponent = 1, nComponents - 1
  IF( ( NonSphericalMixingRule == 1 .OR. NonSphericalMixingRule == 2 ) .AND. ( cGeometry(cComponent) /= &
  &     cGeometry(cComponent+1) ) ) THEN
    WRITE( *, "(2G0)" ) "Error: Molecular geometries MUST MATCH when non-spherical mixing rules are 1 (aspect ratio rule) ", &
    &                   "or 2 (angle average rule). Exiting..."
    CALL EXIT(  )
  END IF
END DO
READ( 10, * ) Dummy, ZhangCorrectionLogical
READ( 10, * ) Dummy, HigherOrderTPTLogical
READ( 10, * ) Dummy, nHigherOrder
IF( HigherOrderTPTLogical .AND. ( nHigherOrder <= 2 ) ) THEN
  WRITE( *, "(G0)" ) "Error: the number of higher order terms must be greater than 2. Exiting..."
  CALL EXIT(  )
END IF
READ( 10, * ) Dummy, EffPFractionMixingRule
IF( EffPFractionMixingRule < 1 .OR. EffPFractionMixingRule > 2 ) THEN
  WRITE( *, "(2G0)" ) "Error: the effective packing fraction mixing rule must be 1 (reduced density 3 mixing rule) or 2 ", &
  &                   "(one-fluid van der Waals mixing rule). Exiting..."
  CALL EXIT(  )
END IF
READ( 10, * ) Dummy, SpecificHeatReference
CALL ToUpper( SpecificHeatReference, LEN_TRIM( SpecificHeatReference ), SpecificHeatReference )
IF( TRIM( SpecificHeatReference ) /= "TRC" .AND. SpecificHeatReference /= "NIST" ) THEN
  WRITE( *, "(2G0)" ) "Error: the specific heat reference must be 'TRC' or 'NIST'. Exiting..."
  CALL EXIT(  )
END IF
READ( 10, * ) Dummy, PotentialType
IF( PotentialType < 1 .OR. PotentialType > 4 ) THEN
  WRITE( *, "(2G0)" ) "Error: the potential type must be 1 (Square-Well), 2 (Sutherland), 3 (Yukawa), or 4 (Convex Square-Well)", &
  &                   ". Exiting..."
  CALL EXIT(  )
END IF
PotentialTypeLogical = .FALSE.
IF( PotentialType == 1 ) THEN
  PotentialTypeLogical(1) = .TRUE.
ELSE IF( PotentialType == 2 ) THEN
  PotentialTypeLogical(2) = .TRUE.
ELSE IF( PotentialType == 3 ) THEN
  PotentialTypeLogical(3) = .TRUE.
ELSE IF( PotentialType == 4 ) THEN
  PotentialTypeLogical(4) = .TRUE.
END IF
READ( 10, * ) Dummy, PYHCBCorrectionLogical
READ( 10, * ) Dummy, UseA1ForA2Logical
READ( 10, * ) Dummy, ReferenceBoublikLogical
READ( 10, * ) Dummy, UseA1AllGeometries
READ( 10, * ) Dummy, UseA2AllGeometries
CLOSE( 10 )

RETURN

END SUBROUTINE Read_Global_Specification

! ************************************************************************************************ !
!                 This subroutine reads the input file and stores the information                  !
! ************************************************************************************************ !
SUBROUTINE Read_Mixture_Specification( mFraction )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction ! Molar fraction

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 01 ) :: Dummy ! Dummy

! Read specifications
OPEN( Unit= 10, File= "specs_mixture.ini", Action= "READ" )
READ( 10, * ) Dummy, mFraction
READ( 10, * ) Dummy, ijWellDepthCorrection
READ( 10, * ) Dummy, ijPotentialRangeCorrection
CLOSE( 10 )

RETURN

END SUBROUTINE Read_Mixture_Specification

! ************************************************************************************************ !
!               This subroutine creates progress bar for the thermodynamic processes               !
! ************************************************************************************************ !
SUBROUTINE Progress_Bar_Path_01( Process, FirstVariable, SecondVariable )

! Use kind Real64 and Int64, and Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: Process ! Process type

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: FirstVariable  ! First variable
REAL( Kind= Real64 ) :: SecondVariable ! Second variable
REAL( Kind= Real64 ) :: Pressure       ! Pressure
REAL( Kind= Real64 ) :: Temperature    ! Temperature
REAL( Kind= Real64 ) :: Density        ! Density

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 66 ) :: ProgressBar ! Progress bar

! Progress bar (FORMAT)
IF( Process == 1 ) THEN ! Isothermal process
  ! Initialization
  Temperature = FirstVariable
  Pressure = SecondVariable
  ! Conversion from Pa to MPa
  Pressure = Pressure / 1.D6
  ! Print progress bar
  ProgressBar = "Current Temperature: ???????? K | Current Pressure: ??????? MPa"
  WRITE( Unit= ProgressBar(22:29), Fmt= "(F8.3)" ) Temperature
  WRITE( Unit= ProgressBar(53:59), Fmt= "(F7.2)" ) Pressure
  ! Print progress bar
  WRITE( Unit= Output_Unit, Fmt= "(A1,A63)", Advance= "NO" ) CHAR(13), ProgressBar
ELSE IF( Process == 2 ) THEN ! Isobaric process
  ! Initialization
  Pressure = FirstVariable
  Temperature = SecondVariable
  ! Conversion from Pa to MPa
  Pressure = Pressure / 1.D6
  ! Print progress bar
  ProgressBar = "Current Pressure: ??????? MPa | Current Temperature: ???????? K"
  WRITE( Unit= ProgressBar(19:25), Fmt= "(F7.2)" ) Pressure
  WRITE( Unit= ProgressBar(54:61), Fmt= "(F8.3)" ) Temperature
  ! Print progress bar
  WRITE( Unit= Output_Unit, Fmt= "(A1,A63)", Advance= "NO" ) CHAR(13), ProgressBar
ELSE IF( Process == 3 ) THEN ! Isochoric process
  ! Initialization
  Density = FirstVariable
  Temperature = SecondVariable
  ! Print progress bar
  ProgressBar = "Current Density: ???????? mol/m3 | Current Temperature: ???????? K"
  WRITE( Unit= ProgressBar(18:25), Fmt= "(F8.1)" ) Density
  WRITE( Unit= ProgressBar(57:64), Fmt= "(F8.3)" ) Temperature
  ! Print progress bar
  WRITE( Unit= Output_Unit, Fmt= "(A1,A66)", Advance= "NO" ) CHAR(13), ProgressBar
END IF

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE Progress_Bar_Path_01

! ************************************************************************************************ !
!               This subroutine creates progress bar for the thermodynamic processes               !
! ************************************************************************************************ !
SUBROUTINE Progress_Bar_Path_02( FirstVariable )

! Use kind Real64 and Int64, and Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: FirstVariable ! First variable
REAL( Kind= Real64 ) :: Temperature   ! Temperature

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 34 ) :: ProgressBar ! Progress bar

! Initialization
Temperature = FirstVariable
! Print progress bar
ProgressBar = "Parsing temperature of: ???????? K"
WRITE( Unit= ProgressBar(25:32), Fmt= "(F8.3)" ) Temperature
! Print progress bar
WRITE( Unit= Output_Unit, Fmt= "(A1,A34)", Advance= "NO" ) CHAR(13), ProgressBar

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE Progress_Bar_Path_02

! ************************************************************************************************ !
!               This subroutine creates progress bar for the thermodynamic processes               !
! ************************************************************************************************ !
SUBROUTINE Progress_Bar_Path_03( FirstVariable, SecondVariable, ThirdVariable, SimplexType, FourthVariable, FifthVariable )

! Use 'nSimplexCycles' and kind Real64 and Int64, and Output_Unit
USE GlobalVariables, ONLY: nSimplexCycles
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: FirstVariable  ! First variable
INTEGER( Kind= Int64 ) :: SecondVariable ! Second variable
INTEGER( Kind= Int64 ) :: ThirdVariable  ! Third variable
INTEGER( Kind= Int64 ) :: iPoint         ! Counter
INTEGER( Kind= Int64 ) :: iGuess         ! Counter
INTEGER( Kind= Int64 ) :: nPoints        ! Total number of points
INTEGER( Kind= Int64 ) :: SimplexType    ! Simplex type (-1: Initial Guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction; 5: Best Objective Function)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: CurrentPercentage ! Current progress
REAL( Kind= Real64 ) :: ObjectiveFunction ! Objective function
REAL( Kind= Real64 ) :: Convergence       ! Convergence criterion
REAL( Kind= Real64 ) :: FourthVariable    ! Fourth variable
REAL( Kind= Real64 ) :: FifthVariable     ! Fifth variable

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 147 ) :: ProgressBar ! Progress bar

! Initialization
IF( SimplexType == -1 ) THEN
  iPoint  = FirstVariable
  nPoints = SecondVariable
  iGuess  = ThirdVariable
  ! Progress
  CurrentPercentage = DBLE( iPoint ) / DBLE( nPoints ) * 100.D0
  ! Print progress bar
  IF( CurrentPercentage >= 0.D0 .AND. CurrentPercentage < 10.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??. "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A31)", Advance= "NO" ) CHAR(13), ProgressBar(1:31)
  ELSE IF( CurrentPercentage >= 10.D0 .AND. CurrentPercentage < 20.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??.. "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A32)", Advance= "NO" ) CHAR(13), ProgressBar(1:32)
  ELSE IF( CurrentPercentage >= 20.D0 .AND. CurrentPercentage < 30.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A33)", Advance= "NO" ) CHAR(13), ProgressBar(1:33)
  ELSE IF( CurrentPercentage >= 30.D0 .AND. CurrentPercentage < 40.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??.... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A34)", Advance= "NO" ) CHAR(13), ProgressBar(1:34)
  ELSE IF( CurrentPercentage >= 40.D0 .AND. CurrentPercentage < 50.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??..... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A35)", Advance= "NO" ) CHAR(13), ProgressBar(1:35)
  ELSE IF( CurrentPercentage >= 50.D0 .AND. CurrentPercentage < 60.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??...... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A36)", Advance= "NO" ) CHAR(13), ProgressBar(1:36)
  ELSE IF( CurrentPercentage >= 60.D0 .AND. CurrentPercentage < 70.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??....... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A37)", Advance= "NO" ) CHAR(13), ProgressBar(1:37)
  ELSE IF( CurrentPercentage >= 70.D0 .AND. CurrentPercentage < 80.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??........ "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A38)", Advance= "NO" ) CHAR(13), ProgressBar(1:38)
  ELSE IF( CurrentPercentage >= 80.D0 .AND. CurrentPercentage < 90.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??......... "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A39)", Advance= "NO" ) CHAR(13), ProgressBar(1:39)
  ELSE IF( CurrentPercentage >= 90.D0 .AND. CurrentPercentage < 100.D0 ) THEN
    ProgressBar = "Calculating Initial Guess #??.........."
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A39)", Advance= "NO" ) CHAR(13), ProgressBar(1:39)
  ELSE
    ProgressBar = "Calculating Initial Guess #0           "
    WRITE( Unit= ProgressBar(28:29), Fmt= "(I0.2)" ) iGuess
    ! Print progress bar
    WRITE( Unit= Output_Unit, Fmt= "(A1,A39)", Advance= "NO" ) CHAR(13), ProgressBar(1:39)
  END IF
  WRITE( Unit= ProgressBar(39:39), Fmt= "(A1)" ) " "
  ! Print progress bar
  WRITE( Unit= Output_Unit, Fmt= "(A1,A39)", Advance= "NO" ) CHAR(13), ProgressBar(1:39)
ELSE
  iPoint  = FirstVariable
  nPoints = SecondVariable
  ObjectiveFunction = FourthVariable
  Convergence = FifthVariable
  ! Progress
  CurrentPercentage = DBLE( iPoint ) / DBLE( nPoints ) * 100.D0
  ! Print progress bar
  ProgressBar(1:81)   = "Simplex Cycle: #?????? | Simplex Type: ??????????? | Current Progress: ??????% | "
  ProgressBar(82:147) = "Best Objective Function:  ????????% | Convergence (STD): ?????????"
  ! Simplex cycle
  IF( nSimplexCycles > 999999 ) THEN
    WRITE( Unit= ProgressBar(16:22), Fmt= "(A7)" ) ">999999"
  ELSE
    WRITE( Unit= ProgressBar(17:22), Fmt= "(I0.6)" ) nSimplexCycles
  END IF
  ! Simplex type
  IF( SimplexType == 0 ) THEN
    WRITE( Unit= ProgressBar(40:50), Fmt= "(A11)" ) "Midpoint   "
  ELSE IF( SimplexType == 1 ) THEN
    WRITE( Unit= ProgressBar(40:50), Fmt= "(A11)" ) "Reflection "
  ELSE IF( SimplexType == 2 ) THEN
    WRITE( Unit= ProgressBar(40:50), Fmt= "(A11)" ) "Expansion  "
  ELSE IF( SimplexType == 3 ) THEN
    WRITE( Unit= ProgressBar(40:50), Fmt= "(A11)" ) "Contraction"
  ELSE IF( SimplexType == 4 ) THEN
    WRITE( Unit= ProgressBar(40:50), Fmt= "(A11)" ) "Reduction  "
  END IF
  ! Progress
  WRITE( Unit= ProgressBar(72:77), Fmt= "(F6.2)" ) CurrentPercentage
  ! Objective function
  IF( ObjectiveFunction > 999.9999D0 ) THEN
    WRITE( Unit= ProgressBar(108:115), Fmt= "(A8)" ) ">999.999"
  ELSE
    WRITE( Unit= ProgressBar(108:115), Fmt= "(F8.4)" ) ObjectiveFunction
  END IF
  ! Convergence
  IF( Convergence > 999.9999D0 ) THEN
    WRITE( Unit= ProgressBar(139:147), Fmt= "(A9)" ) ">999.9999"
  ELSE IF( Convergence >= 1.D0 ) THEN
    WRITE( Unit= ProgressBar(139:147), Fmt= "(F8.4)" ) Convergence
  ELSE
    WRITE( Unit= ProgressBar(139:147), Fmt= "(E9.3)" ) Convergence
  END IF
  ! Print progress bar
  WRITE( Unit= Output_Unit, Fmt= "(A1,A147)", Advance= "NO" ) CHAR(13), ProgressBar(1:147)
END IF

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE Progress_Bar_Path_03

! ************************************************************************************************ !
!           This subroutine gets the current date and time to organize files and folders           !
! ************************************************************************************************ !
SUBROUTINE Get_Date_and_Time(  )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER, DIMENSION( 8 ) :: DateTime ! Date and time

! Get current date and time
CALL DATE_AND_TIME( VALUES= DateTime )

! Date format (YYYY/MM/DD)
FormatDate = "(I4,2I2.2)"
WRITE( DescriptorDate, FormatDate ) DateTime(1), DateTime(2), DateTime(3)

! Time format (HH:MM:SS)
FormatHour = "(3I2.2)"
WRITE( DescriptorHour, FormatHour ) DateTime(5), DateTime(6), DateTime(7)

RETURN

END SUBROUTINE Get_Date_and_Time

! ************************************************************************************************ !
!               This subroutine sorts an array in ascending order using the bubble sort            !
! ************************************************************************************************ !
SUBROUTINE Bubble_Sort( ArraySize, Array, ArrayAux )

! Uses one module: global variables
USE GlobalVariables, ONLY: Int64, Real64

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iArray, jArray ! Counters
INTEGER( Kind= Int64 ) :: ArraySize      ! Array size

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                                      :: ArrayTemp    ! Temporary array
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

! ************************************************************************************************ !
!               This subroutine sorts an array in ascending order using the bubble sort            !
! ************************************************************************************************ !
SUBROUTINE Bubble_Sort_VLE( ArraySize, Array )

! Uses one module: global variables
USE GlobalVariables, ONLY: Int64, Real64

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iArray, jArray ! Counters
INTEGER( Kind= Int64 ) :: ArraySize      ! Array size

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                         :: ArrayTemp ! Temporary variable
REAL( Kind= Real64 ), DIMENSION( ArraySize ) :: Array     ! Array

! Sort array
DO iArray = 1, ArraySize-1
  DO jArray = 1, ArraySize-iArray
    IF( Array(jArray) > Array(jArray+1) ) THEN
      ArrayTemp        = Array(jArray)
      Array(jArray)    = Array(jArray+1)
      Array(jArray+1)  = ArrayTemp
    END IF
  END DO
END DO

RETURN

END SUBROUTINE Bubble_Sort_VLE