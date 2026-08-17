! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!         This module creates folders and subfolders to organize the simulation results.          !
! Directories are created by executing a shell command via an intrinsic function called 'system'. !
!               Please note that which shell is used to invoke the command line is                !
!                           system-dependent and environment-dependent.                           !
!         See <https://gcc.gnu.org/onlinedocs/gfortran/SYSTEM.html> for more information.         !
!                     The code below is meant for Linux operational systems.                      !
!                   We have not provided an alternative code for Windows users.                   !
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
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE Folders

! Uses one MODULE: global variables
USE GlobalVariables

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!                              INITIALIZATION OF PARENT DIRECTORIES                               !
! *********************************************************************************************** !
SUBROUTINE InitialConfigFolder(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration", Exist= FolderExistLogical(1) )

! Initial configuration folder (holds information on the initial molecular structure)
IF( .NOT. FolderExistLogical(1) ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
! The initial molecular structure at 'OVITO' subfolder is properly formatted to be analyzed by that software.
INQUIRE( File= "Initial_Configuration/OVITO/", Exist= SubfolderExistLogical(1) )

! Initial configuration subfolder
IF( .NOT. SubfolderExistLogical(1) ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/" )
END IF

! Date format (YYYY/MM/DD)
DateFormat = "(I4,2I2.2)"
! Date descriptor
WRITE( DateDescriptor, DateFormat ) DateTime(1), DateTime(2), DateTime(3)
! Time format (HH:MM:SS)
HourFormat = "(3I2.2)"
! Hour descriptor
WRITE( HourDescriptor, HourFormat ) DateTime(5), DateTime(6), DateTime(7)

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration/OVITO/"//TRIM( DateDescriptor )//"/", Exist= DateFolderExistLogical(1) )

! Date subfolder
IF( .NOT. DateFolderExistLogical(1) ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/"//TRIM( DateDescriptor )//"/" )
END IF

RETURN

END SUBROUTINE InitialConfigFolder

! *********************************************************************************************** !
!                              INITIALIZATION OF PARENT DIRECTORIES                               !
! *********************************************************************************************** !
SUBROUTINE MainFolders(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Trajectories", Exist= FolderExistLogical(2) )
INQUIRE( File= "Ratio", Exist= FolderExistLogical(3) )
INQUIRE( File= "Results", Exist= FolderExistLogical(4) )
INQUIRE( File= "Floppy-box", Exist= FolderExistLogical(5) )

! Trajectory folder (holds information on orientation and position of particles)
IF( .NOT. FolderExistLogical(2) ) THEN
  CALL SYSTEM( "mkdir Trajectories" )
END IF

! Ratio folder (holds information on the equilibration cycles, like maximum displacement adjustment)
IF( .NOT. FolderExistLogical(3) ) THEN
  CALL SYSTEM( "mkdir Ratio" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Translation/", Exist= SubfolderExistLogical(2) )
INQUIRE( File= "Ratio/Rotation/", Exist= SubfolderExistLogical(3) )

! Ratio subfolders
IF( .NOT. SubfolderExistLogical(2) ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/" )
END IF
IF( .NOT. SubfolderExistLogical(3) ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/" )
END IF

! Results folder (holds information on the packing fraction of the system)
IF( .NOT. FolderExistLogical(4) ) THEN
  CALL SYSTEM( "mkdir Results" )
END IF

! Floppy box folder (holds information on properties of the floppy box)
IF( .NOT. FolderExistLogical(5) ) THEN
  CALL SYSTEM( "mkdir Floppy-box" )
END IF

RETURN

END SUBROUTINE MainFolders

! *********************************************************************************************** !
!                                INITIALIZATION OF DATE SUBFOLDERS                                !
! *********************************************************************************************** !
SUBROUTINE DateFolders(  )

IMPLICIT NONE

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Trajectories/"//TRIM( DateDescriptor )//"/", Exist= DateFolderExistLogical(1) )
INQUIRE( File= "Ratio/Translation/"//TRIM( DateDescriptor )//"/", Exist= DateFolderExistLogical(2) )
INQUIRE( File= "Ratio/Rotation/"//TRIM( DateDescriptor )//"/", Exist= DateFolderExistLogical(3) )
INQUIRE( File= "Results/"//TRIM( DateDescriptor )//"/", Exist= DateFolderExistLogical(4) )

! Date subfolders
IF( .NOT. DateFolderExistLogical(1) ) THEN
  CALL SYSTEM( "mkdir Trajectories/"//TRIM( DateDescriptor )//"/" )
END IF
IF( .NOT. DateFolderExistLogical(2) ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/"//TRIM( DateDescriptor )//"/" )
END IF
IF( .NOT. DateFolderExistLogical(3) ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/"//TRIM( DateDescriptor )//"/" )
END IF
IF( .NOT. DateFolderExistLogical(4) ) THEN
  CALL SYSTEM( "mkdir Results/"//TRIM( DateDescriptor )//"/" )
END IF

RETURN

END SUBROUTINE DateFolders

END MODULE Folders