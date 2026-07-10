! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!         This module creates folders and subfolders to organize the simulation results.          !
! Directories are created by executing a shell command via an intrinsic function called 'system'. !
!               Please note that which shell is used to invoke the command line is                !
!                           system-dependent and environment-dependent.                           !
!         See <https://gcc.gnu.org/onlinedocs/gfortran/SYSTEM.html> for more information.         !
!                  The code below is meant for Ubuntu 20.04 operational systems.                  !
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
!                                       February 15th, 2024                                       !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE Folders

! Uses one MODULE: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FolderExist     ! Checks whether folder exists or not
LOGICAL :: DateFolderExist ! Checks whether date folders exist or not
LOGICAL :: SubfolderExist  ! Checks whether subfolder exists or not
LOGICAL :: AdvanceLine     ! Checks whether subfolder exists or not

CONTAINS

! *********************************************************************************************** !
!                             Initialization of parental directories                              !
! *********************************************************************************************** !
SUBROUTINE InitialConfigurationFolder(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration", Exist= FolderExist )

! Initial configuration folder (holds information on the initial molecular structure)
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating initial configuration folder..."
  CALL System( "mkdir Initial_Configuration" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration/OVITO/", Exist= SubfolderExist )

! Initial configuration subfolder
IF( .NOT. SubfolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating initial configuration subfolder..."
  CALL System( "mkdir Initial_Configuration/OVITO/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )

! Date subfolder
IF( .NOT. DateFolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolder..."
  CALL System( "mkdir Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE InitialConfigurationFolder

! *********************************************************************************************** !
!                             Initialization of parental directories                              !
! *********************************************************************************************** !
SUBROUTINE InquireFolders(  )

IMPLICIT NONE

! Initialization
AdvanceLine = .FALSE.

! Trajectory folder (holds information on orientation and position of particles)
INQUIRE( File= "Trajectories/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating trajectory folder..."
  CALL System( "mkdir Trajectories/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Ratio folder (holds information on the equilibration cycles)
INQUIRE( File= "Ratio/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating ratio folder..."
  CALL System( "mkdir Ratio/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Box folder (holds information on the box properties)
INQUIRE( File= "Box/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating box folder..."
  CALL System( "mkdir Box/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Order parameter folder (holds information on the nematic order parameter)
INQUIRE( File= "Order_Parameter/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating order parameter folder..."
  CALL System( "mkdir Order_Parameter/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Results folder (holds information on the packing fraction and the simulation box)
INQUIRE( File= "Results/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating results folder..."
  CALL System( "mkdir Results/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Inquires whether a subfolder exists or not
INQUIRE( File= "Ratio/Translation/", Exist= SubfolderExist )
IF( .NOT. SubfolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating ratio subfolders..."
  AdvanceLine = .TRUE.
ELSE
  INQUIRE( File= "Ratio/Rotation/", Exist= SubfolderExist )
  IF( .NOT. SubfolderExist ) THEN
    WRITE( *, "(G0)", Advance= "NO" ) "Creating ratio subfolders..."
    AdvanceLine = .TRUE.
  ELSE
    INQUIRE( File= "Ratio/Volume/", Exist= SubfolderExist )
    IF( .NOT. SubfolderExist ) THEN
      WRITE( *, "(G0)", Advance= "NO" ) "Creating ratio subfolders..."
      AdvanceLine = .TRUE.      
    END IF
  END IF
END IF

! Translational ratio subfolder
INQUIRE( File= "Ratio/Translation/", Exist= SubfolderExist )
IF( .NOT. SubfolderExist ) THEN
  CALL System( "mkdir Ratio/Translation/" )
END IF

! Rotational ratio subfolder
INQUIRE( File= "Ratio/Rotation/", Exist= SubfolderExist )
IF( .NOT. SubfolderExist ) THEN
  CALL System( "mkdir Ratio/Rotation/" )
END IF

! Volumetric ratio subfolder
INQUIRE( File= "Ratio/Volume/", Exist= SubfolderExist )
IF( .NOT. SubfolderExist ) THEN
  CALL System( "mkdir Ratio/Volume/" )
END IF

! Skip line
IF( AdvanceLine ) THEN
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE InquireFolders

! *********************************************************************************************** !
!                                Initialization of date subfolders                                !
! *********************************************************************************************** !
SUBROUTINE DateFolders(  )

IMPLICIT NONE

! Initialization
AdvanceLine = .FALSE.

! Inquires whether the date subfolder exists or not
INQUIRE( File= "Trajectories/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
  AdvanceLine = .TRUE.
ELSE
  INQUIRE( File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
  IF( .NOT. DateFolderExist ) THEN
    WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
    AdvanceLine = .TRUE.
  ELSE
    INQUIRE( File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
    IF( .NOT. DateFolderExist ) THEN
      WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
      AdvanceLine = .TRUE.
    ELSE
      INQUIRE( File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
      IF( .NOT. DateFolderExist ) THEN
        WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
        AdvanceLine = .TRUE.
      ELSE
        INQUIRE( File= "Box/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
        IF( .NOT. DateFolderExist ) THEN
          WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
          AdvanceLine = .TRUE.
        ELSE
          INQUIRE( File= "Order_Parameter/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
          IF( .NOT. DateFolderExist ) THEN
            WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
            AdvanceLine = .TRUE.
          ELSE
            INQUIRE( File= "Results/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
            IF( .NOT. DateFolderExist ) THEN
              WRITE( *, "(G0)", Advance= "NO" ) "Creating date subfolders..."
              AdvanceLine = .TRUE.
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
END IF

! Trajectory date subfolder
INQUIRE( File= "Trajectories/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Trajectories/"//TRIM( DescriptorDate )//"/" )
END IF

! Ratios date subfolder
INQUIRE( File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Ratio/Translation/"//TRIM( DescriptorDate )//"/" )
END IF
INQUIRE( File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Ratio/Rotation/"//TRIM( DescriptorDate )//"/" )
END IF
INQUIRE( File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Ratio/Volume/"//TRIM( DescriptorDate )//"/" )
END IF

! Box date subfolder
INQUIRE( File= "Box/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Box/"//TRIM( DescriptorDate )//"/" )
END IF

! Order parameter date subfolder
INQUIRE( File= "Order_Parameter/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Order_Parameter/"//TRIM( DescriptorDate )//"/" )
END IF

! Results date subfolder
INQUIRE( File= "Results/"//TRIM( DescriptorDate )//"/", Exist= DateFolderExist )
IF( .NOT. DateFolderExist ) THEN
  CALL System( "mkdir Results/"//TRIM( DescriptorDate )//"/" )
END IF

! Skip line
IF( AdvanceLine ) THEN
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE DateFolders

! *********************************************************************************************** !
!                               Initialization of the backup folder                               !
! *********************************************************************************************** !
SUBROUTINE BackupFolder(  )

IMPLICIT NONE

! Inquires whether the backup folder exists or not
INQUIRE( File= "Backup/", Exist= FolderExist )
IF( .NOT. FolderExist ) THEN
  WRITE( *, "(G0)", Advance= "NO" ) "Creating backup folder..."
  CALL System( "mkdir Backup/" )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE BackupFolder

END MODULE Folders