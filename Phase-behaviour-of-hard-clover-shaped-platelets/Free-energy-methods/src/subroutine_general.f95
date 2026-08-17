! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!              This code contains the general subroutines used in the main program.               !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                        August 13th, 2026                                        !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
!                             ---------------------------------------                             !
!                     J. Graaf, L. Filion, M. Marechal, R. Roij, M. Dijkstra                      !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                  This subroutine generates a progress bar for the simulation.                   !
! *********************************************************************************************** !
SUBROUTINE ProgressBar( CurrentCycle, TotalCycles, EnsembleType )

! Uses two modules
USE GlobalVariables

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: CurrentCycle ! Counter
INTEGER( Kind= Int64 ) :: TotalCycles  ! Total/Maximum number of cycles

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 21 ) :: ProgressBarMC ! Progress bar
CHARACTER( LEN= 03 ) :: EnsembleType  ! Ensemble type

! Progress bar (FORMAT)
IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 10.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ???%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 100.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ????%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 1000.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ?????%"
END IF

! Progress bar (replace character positions)
IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 10.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:18), FMT= "(F3.1)" ) ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 - 0.05D0
  IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 - 0.05D0 < 0.D0 ) THEN
    WRITE( Unit= ProgressBarMC(16:18), FMT= "(F3.1)" ) 0.D0
  END IF
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 100.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:19), FMT= "(F4.1)" ) ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 - 0.05D0
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 1000.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16:20), FMT= "(F5.1)" ) ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0
END IF

! Print progress bar
IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 10.D0 ) THEN
  WRITE( Unit= Output_Unit, FMT= "(A1,A19)" , ADVANCE= "NO" ) CHAR(13), ProgressBarMC
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 100.D0 ) THEN
  WRITE( Unit= Output_Unit, FMT= "(A1,A20)" , ADVANCE= "NO" ) CHAR(13), ProgressBarMC
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 1000.D0 ) THEN
  WRITE( Unit= Output_Unit, FMT= "(A1,A21)" , ADVANCE= "NO" ) CHAR(13), ProgressBarMC
END IF

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBar

! *********************************************************************************************** !
!                 This function converts any string into uppercase (from A to Z)                  !
! *********************************************************************************************** !
SUBROUTINE ToUpper( StringInput, StringLength, StringOutput )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: StringLength ! String length
INTEGER :: iCharacter   ! ASCII character code
INTEGER :: iString      ! Counter

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
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