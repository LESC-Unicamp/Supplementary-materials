MODULE Folders

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!                              INITIALIZATION OF PARENT DIRECTORIES                               !
! *********************************************************************************************** !
SUBROUTINE Initialize_Folders(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration", Exist= FolderExist )

! Initial configuration folder (holds information on the initial molecular structure)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
! The initial molecular structure at 'OVITO' subfolder is properly formatted to be analyzed by that software.
INQUIRE( File= "Initial_Configuration/OVITO/", Exist= FolderExist )

! Initial configuration subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Backup/", Exist= FolderExist )

! Backup folder (holds information on the backup of the system)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Backup" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Trajectories", Exist= FolderExist )

! Trajectory folder (holds information on orientation and position of particles)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Trajectories" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio", Exist= FolderExist )

! Ratio folder (holds information on the equilibration cycles, like maximum displacement adjustment)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Translation/", Exist= FolderExist )

! Translation subfolders
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Rotation/", Exist= FolderExist )

! Rotation subfolders
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Volume/", Exist= FolderExist )

! Volume subfolders
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Volume/" )
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Box/", Exist= FolderExist )

! Box subfolders
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Box/" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Results", Exist= FolderExist )

! Results folder (holds information on the packing fraction of the system)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Results" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Floppy-box", Exist= FolderExist )

! Floppy box folder (holds information on properties of the floppy box)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Floppy-box" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Curve", Exist= FolderExist )

! Curve folder (holds information on the pressure-density curve)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Curve" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Configurations", Exist= FolderExist )

! Configurations folder (holds information on the configurations of previous simulated systems)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Configurations" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Potential", Exist= FolderExist )

! Potential folder (holds information on the potential of the system)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Potential" )
END IF

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Stack_Rotation", Exist= FolderExist )

! Stack rotation folder (holds information on the rotation of the stack)
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Stack_Rotation" )
END IF

! Date format (YYYY/MM/DD)
FormatDate = "(I4,2I2.2)"
! Date descriptor
WRITE( DescriptorDate, FormatDate ) DateTime(1), DateTime(2), DateTime(3)
! Time format (HH:MM:SS)
FormatHour = "(3I2.2)"
! Hour descriptor
WRITE( DescriptorHour, FormatHour ) DateTime(5), DateTime(6), DateTime(7)

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Trajectories/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Trajectories/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Translation/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Rotation/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Ratio/Volume/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Ratio/Volume/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Box/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Box/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Results/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Results/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Potential/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Potential/"//TRIM(DescriptorDate)//"/" )
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( File= "Stack_Rotation/"//TRIM(DescriptorDate)//"/", Exist= FolderExist )

! Date subfolder
IF( .NOT. FolderExist ) THEN
  CALL SYSTEM( "mkdir Stack_Rotation/"//TRIM(DescriptorDate)//"/" )
END IF

RETURN

END SUBROUTINE Initialize_Folders

END MODULE Folders