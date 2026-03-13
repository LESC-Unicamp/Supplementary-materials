! ************************************************************************************* !
!                                     PROGRESS BAR                                      !
! ************************************************************************************* !
!          This subroutine generates a progress bar for the simplex algorithm.          !
! ************************************************************************************* !
SUBROUTINE Progress_Bar( CurrentCycle, TotalCycles, EnsembleType )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: CurrentCycle ! Counter (cycles)
INTEGER( Kind= Int64 ) :: TotalCycles  ! Total/Maximum number of cycles

! CHARACTER STRINGS
CHARACTER( Len= 21 ) :: ProgressBarMC ! Progress bar
CHARACTER( Len= 03 ) :: EnsembleType  ! Ensemble type

! Progress bar (FORMAT)
IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 10.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ???%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 100.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ????%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 1000.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ?????%"
END IF

! Progress bar (replace character positions)
IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 10.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:18), Fmt= "(F3.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 - 0.05D0
  IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 - 0.05D0 < 0.D0 ) THEN
    WRITE( Unit= ProgressBarMC(16:18), Fmt= "(F3.1)" ) 0.D0
  END IF
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 100.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:19), Fmt= "(F4.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 - 0.05D0
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 1000.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:20), Fmt= "(F5.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0
END IF

! Print progress bar
IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 10.D0 ) THEN
  WRITE( Unit= Output_Unit, Fmt= "(A1,A19)" , Advance= "NO" ) CHAR(13), ProgressBarMC
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 100.D0 ) THEN
  WRITE( Unit= Output_Unit, Fmt= "(A1,A20)" , Advance= "NO" ) CHAR(13), ProgressBarMC
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 1000.D0 ) THEN
  WRITE( Unit= Output_Unit, Fmt= "(A1,A21)" , Advance= "NO" ) CHAR(13), ProgressBarMC
END IF

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE Progress_Bar

! *********************************************************************************************** !
!                 This function converts any string into uppercase (from A to Z)                  !
! *********************************************************************************************** !
SUBROUTINE ToUpper( StringInput, StringLength, StringOutput )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER :: StringLength ! String length
INTEGER :: iCharacter   ! ASCII character code
INTEGER :: iString      ! Counter

! CHARACTER STRINGS
CHARACTER( Len= StringLength ) :: StringInput  ! Length of input string
CHARACTER( Len= StringLength ) :: StringOutput ! Length of output string

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