! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!           This module allows the user to choose the initial molecular configuration.            !
!  The only configurations available at the moment are: floppy box (FB) and liquid crystal (LC).  !
!        Molecules will be then allocated in accordance to the selected crystal structure.        !
!     This module also writes out a file containing all particles' positions and quaternions.     !
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
! Main References:                 M. P. Allen, D. J. Tildesley                                   !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             doi: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE InitialConfiguration

! Uses three modules: global variables, quaternion operations, and vector operations
USE GlobalVariables
USE QuaternionOperations
USE VectorOperations

IMPLICIT NONE

! *********************************************************************************************** !
!                                      INITIAL CONFIGURATION                                      !
! *********************************************************************************************** !
CONTAINS

! *********************************************************************************************** !
!          This subroutine allows the user to choose the initial molecular configuration          !
! *********************************************************************************************** !
SUBROUTINE ConfigurationSelection(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = FB                                                                                       !
!  (2) = LC                                                                                       !
! *********************************************************************************************** !
ConfigSelectionLogical(:) = .FALSE.

! Molecular properties
OPEN( Unit= 100, File= "ini_montecarlo.ini", Action= "READ" )
! Molecular structure
READ( 100, * ) DummyText, ConfigName
CALL ToUpper( ConfigName, LEN_TRIM( ConfigName ), ConfigName )
IF( ConfigName == "FB" ) THEN ! Floppy box
  ConfigSelectionLogical(1) = .TRUE.
ELSE IF( ConfigName == "LC" ) THEN ! Liquid crystal
  ConfigSelectionLogical(2) = .TRUE.
END IF
CLOSE( 100 )

RETURN

END SUBROUTINE ConfigurationSelection

! *********************************************************************************************** !
!         This subroutine allows the user to choose the surface geometry of the molecules         !
! *********************************************************************************************** !
SUBROUTINE GeometrySelection(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: FrameLEFT  ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT ! Box frame dimension

! *********************************************************************************************** !
!  FOUR-LEAF CLOVER GEOMETRIES                                                                    !
!                                                                                                 !
!  The geometry of the cross sections consists of four circles and are constructed as follows:    !
!                                                                                                 !
!  (1) = RULES: A circle intersects the centers of the two nearest neighboring circles            !
!               A circle is equidistant from the two nearest neighboring circles                  !
!                                                                                                 !
!  (2) = RULES: All four circles intersect the center of mass of the composed geometry            !
!               A circle is equidistant from the two nearest neighboring circles                  !
! *********************************************************************************************** !
GeometrySelectionLogical(:) = .FALSE.

IF( ConfigSelectionLogical(1) ) THEN
  ! Floppy-box folder
  INQUIRE( File= "Floppy-box", Exist= FolderExistLogical(1) )
  ! Trajectory folder (holds information on orientation and position of particles)
  IF( .NOT. FolderExistLogical(1) ) THEN
    CALL SYSTEM( "mkdir Floppy-box" )
  END IF
  ! Floppy-box file
  INQUIRE( File= "Floppy-box/fbox.dat", Exist= FileExistLogical )
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  FrameLEFT  = (35 - FLOOR( REAL( 15 ) / 2.D0 ) )
  FrameRIGHT = (35 - CEILING( REAL( 15 ) / 2.D0 ) )
  WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"FLOPPY-BOX DATA"//REPEAT( " ", FrameRIGHT )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  IF( .NOT. FileExistLogical ) THEN
    WRITE( *, "(G0)" ) "Floppy-box data file 'Floppy-box/fbox.dat' not found! See 'README.md' for instructions. Exiting..."
    WRITE( *, "(G0)") " "
    STOP
  ELSE IF( FileExistLogical ) THEN
    WRITE( *, "(G0)" ) "Floppy-box data file found! Resuming..."
    WRITE( *, "(G0)" ) " "
  END IF
  ! Molecular properties
  OPEN( Unit= 10, File= "Floppy-box/fbox.dat", Action= "READ" )
  ! Cross-sectional arrangement
  READ( 10, * ) DummyText, cArrangement
  GeometrySelectionLogical(:) = .FALSE.
  ! Arrangement (1)
  IF( cArrangement == 1 ) THEN
    GeometrySelectionLogical(1) = .TRUE.
  ! Arrangement (2)
  ELSE IF( cArrangement == 2 ) THEN
    GeometrySelectionLogical(2) = .TRUE.
  END IF
  ! Cylindrical diameter
  READ( 10, * ) DummyText, cDiameter ! Å
  ! Cylindrical length
  READ( 10, * ) DummyText, cLength ! Å
  ! Aspect ratio of cylinders
  cAspectRatio = cLength / cDiameter
  CLOSE( 10 )
ELSE IF( ConfigSelectionLogical(2) ) THEN
  ! Liquid-crystal folder
  INQUIRE( File= "Liquid-crystal", Exist= FolderExistLogical(1) )
  ! Trajectory folder (holds information on orientation and position of particles)
  IF( .NOT. FolderExistLogical(1) ) THEN
    CALL SYSTEM( "mkdir Liquid-crystal" )
  END IF
  ! Liquid-crystal file
  INQUIRE( File= "Liquid-crystal/lc.dat", Exist= FileExistLogical )
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  FrameLEFT  = (35 - FLOOR( REAL( 15 ) / 2.D0 ) )
  FrameRIGHT = (35 - CEILING( REAL( 15 ) / 2.D0 ) )
  WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"LIQUID-CRYSTAL DATA"//REPEAT( " ", FrameRIGHT )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  IF( .NOT. FileExistLogical ) THEN
    WRITE( *, "(G0)" ) "Liquid-crystal data file 'Liquid-crystal/lc.dat' not found! See 'README.md' for instructions. Exiting..."
    WRITE( *, "(G0)") " "
    STOP
  ELSE IF( FileExistLogical ) THEN
    WRITE( *, "(G0)" ) "Liquid-crystal data file found! Resuming..."
    WRITE( *, "(G0)" ) " "
  END IF
  ! Liquid-crystal data
  OPEN( Unit= 10, File= "Liquid-crystal/lc.dat", Action= "READ" )
  ! Cross-sectional arrangement
  READ( 10, * ) DummyText, cArrangement
  GeometrySelectionLogical(:) = .FALSE.
  ! Arrangement (1)
  IF( cArrangement == 1 ) THEN
    GeometrySelectionLogical(1) = .TRUE.
  ! Arrangement (2)
  ELSE IF( cArrangement == 2 ) THEN
    GeometrySelectionLogical(2) = .TRUE.
  END IF
  ! Cylindrical diameter
  READ( 10, * ) DummyText, cDiameter ! Å
  ! Cylindrical length
  READ( 10, * ) DummyText, cLength ! Å
  ! Aspect ratio of cylinders
  cAspectRatio = cLength / cDiameter
  CLOSE( 10 )
END IF

RETURN

END SUBROUTINE GeometrySelection

! *********************************************************************************************** !
!               Subroutine to calculate the surface area of the non-convex geometry               !
!       The non-convex geometry is composed of four identical cylinders, such that [1] the        !
! circumference of a cylinder passes through the centers of the two nearest cylinders or [2] the  !
!           circumference of all cylinders passes through the molecular center of mass            !
!                                  (common intersection point).                                   !
!   In the case [1], the distance between the centers of the circumferences of the two nearest    !
!     cylinders is the radius of one cylinder, and the width of the non-convex body is 1.5D,      !
!                            where D is the diameter of the cylinder.                             !
!   In the case [2], the distance between the centers of the circumferences of the two nearest    !
!      cylinders is the hypotenuse (hyp) of an isosceles triangle whose side is equal to the      !
!            radius of one cylinder, and the width of the non-convex body is D + hyp.             !
! *********************************************************************************************** !
SUBROUTINE CrossSectionalArea( SurfaceArea )

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: aCircle, aIntersection2Circles, aIntersection4Circles ! Area of the circles and their intersections
REAL( Kind= Real64 ) :: SurfaceArea                                           ! Surface area

! GEOMETRY[1]
IF( GeometrySelectionLogical(1) ) THEN
  ! Area of circles
  aCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  aIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( 2.D0 * cPi / 3.D0 ) - ( DSQRT( 3.D0 ) / 2.D0 ) )
  ! Area of the intersection between four circles
  aIntersection4Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( cPi / 3.D0 ) - ( DSQRT( 3.D0 ) ) + 1.D0 )
  ! Surface area
  SurfaceArea = (4.D0 * aCircle ) - ( 4.D0 * aIntersection2Circles ) + ( aIntersection4Circles )
! GEOMETRY[2]
ELSE IF( GeometrySelectionLogical(2) ) THEN
  ! Area of circles
  aCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  aIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( cPi / 2.D0 - 1.D0 )
  ! Surface area
  SurfaceArea = 4.D0 * aCircle - 4.D0 * aIntersection2Circles
END IF

RETURN

END SUBROUTINE CrossSectionalArea

! *********************************************************************************************** !
!            This subroutine allocates particles according to the FB configuration                !
! *********************************************************************************************** !
SUBROUTINE ConfigFloppyBox(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCell, jCell, kCell ! Counters (cell)
INTEGER( Kind= Int64 ) :: pParticle, pCounter ! Counters (particle)
INTEGER( Kind= Int64 ) :: pCell, ParticleCell ! Counters (particle in unit cell)
INTEGER( Kind= Int64 ) :: nCells              ! Number of unit cells
INTEGER( Kind= Int64 ) :: FrameLEFT           ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT          ! Box frame dimension


! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: CubicRootCheck         ! Cube root verification
REAL( Kind= Real64 )                    :: CellLengthCorrection   ! Scalar density-dependent factor to correct length of unit cell
REAL( Kind= Real64 )                    :: CellVolume             ! Cell volume
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CellLength             ! Length of unit cell (triclinic)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CellLengthInverse      ! Inverse of length of unit cell (triclinic)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox ! Scaling factor (unit cubic box)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis  ! Relative distance (cylinders)

! Floppy-box data
OPEN( Unit= 10, File= "Floppy-box/fbox.dat", Action= "READ" )
! Skip
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
! Cell length
READ( 10, * ) DummyText, CellLength(1:9)
! Cell volume
READ( 10, * ) DummyText, CellVolume
! Skip
READ( 10, * )
READ( 10, * )
! Number of particles in unit cell
READ( 10, * ) DummyText, pCell
! Allocation
ALLOCATE( pPositionCell(3,pCell), pPositionEq(3,nParticles) )
ALLOCATE( pOrientationCell(3,3,pCell), pOrientationEq(3,3,nParticles) )
ALLOCATE( pQuaternionCell(0:3,pCell), pQuaternionEq(0:3,nParticles) )
! Skip
READ( 10, * )
READ( 10, * )
READ( 10, * )
DO pParticle = 1, pCell
  ! Equilibrium positions and orientations (quaternions) in the unit cell
  READ( 10, * ) DummyText, pPositionCell(1:3,pParticle), pQuaternionCell(0:3,pParticle)
  ! Orientation along x-direction
  CALL VectorRotation( xAxis, pQuaternionCell(:,pParticle), pOrientationCell(1,:,pParticle) )
  ! Orientation along y-direction
  CALL VectorRotation( yAxis, pQuaternionCell(:,pParticle), pOrientationCell(2,:,pParticle) )
  ! Orientation along z-direction
  CALL VectorRotation( zAxis, pQuaternionCell(:,pParticle), pOrientationCell(3,:,pParticle) )
END DO
CLOSE( 10 )

! Validation of the number of particles
ParticleValidation: DO
  CubicRootCheck = ( DBLE( nParticles ) / DBLE( pCell ) ) ** ( 1.D0 / 3.D0 )
  IF ( DABS( CubicRootCheck - DNINT( CubicRootCheck ) ) <= 1.D-10 ) THEN
    EXIT ParticleValidation
  ELSE
    WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
    FrameLEFT  = (35 - FLOOR( REAL( 23 ) / 2.D0 ) )
    FrameRIGHT = (35 - CEILING( REAL( 23 ) / 2.D0 ) )
    WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"VALIDATION OF STRUCTURE"//REPEAT( " ", FrameRIGHT )//CH_VS
    WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
    WRITE( *, "(G0)") "Invalid number of particles for the selected floppy-box configuration!"
    WRITE( *, "(5G0)") "The total number of particles (", nParticles,") divided by ", pCell, &
    &                  " (number of particles inside a unit cell) must be a perfect cube root."
    WRITE( *, "(G0)") " "
    WRITE( *, "(G0)") "Exiting..."
    WRITE( *, "(G0)") " "
    STOP
  END IF
END DO ParticleValidation

! Inverse of cell length
CALL InverseMatrixCofactorVec( CellLength, CellLengthInverse, CellVolume )

! Number of unit cells per axis
nCells = NINT( ( DBLE( nParticles ) / DBLE( pCell ) ) ** ( 1.D0 / 3.D0 ) )

! Scalar factor based on the number density of the unit cell
CellLengthCorrection = ( DBLE( pCell ) / ( nDensity * CellVolume ) ) ** ( 1.D0 / 3.D0 )

! New cell length based on the number density of the unit cell
CellLength(:) = CellLengthCorrection * CellLength(:)

! New positions based on the number density of the unit cell
DO pParticle = 1, pCell
  ! Scaled coordinates using old cell length
  CALL MatrixVectorMultiplication( CellLengthInverse, pPositionCell(:,pParticle), ScalingDistanceUnitBox )
  ! New real coordinates using new cell length
  CALL MatrixVectorMultiplication( CellLength, ScalingDistanceUnitBox, pPositionCell(:,pParticle) )
END DO

! New inverse of cell length and volume
CALL InverseMatrixCofactorVec( CellLength, CellLengthInverse, CellVolume )

! Simulation box length
BoxLength(:) = CellLength(:) * DBLE( nCells )

! Inverse of simulation box length
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Position of nonconvex particles (centers of mass)
pCounter = 1
ParticlePlacement: DO
  DO ParticleCell = 1, pCell
    ! Scaled equilibrium position of particle in unit cell
    CALL MatrixVectorMultiplication( CellLengthInverse, pPositionCell(:,ParticleCell), ScalingDistanceUnitBox )
    DO iCell = 1, nCells
      DO jCell = 1, nCells
        DO kCell = 1, nCells
          ! Position of particles in a unit cubic box
          pPosition(1,pCounter) = ScalingDistanceUnitBox(1) + DBLE( iCell - 1 )
          pPosition(2,pCounter) = ScalingDistanceUnitBox(2) + DBLE( jCell - 1 )
          pPosition(3,pCounter) = ScalingDistanceUnitBox(3) + DBLE( kCell - 1 )
          ! Quaternion of particles
          pQuaternion(0,pCounter) = pQuaternionCell(0,ParticleCell)
          pQuaternion(1,pCounter) = pQuaternionCell(1,ParticleCell)
          pQuaternion(2,pCounter) = pQuaternionCell(2,ParticleCell)
          pQuaternion(3,pCounter) = pQuaternionCell(3,ParticleCell)
          ! Equilibrium orientation (quaternion) of particles
          pQuaternionEq(:,pCounter) = pQuaternionCell(:,ParticleCell)
          ! Equilibrium orientation of particles
          pOrientationEq(:,:,pCounter) = pOrientationCell(:,:,ParticleCell)
          ! Increment
          pCounter = pCounter + 1
          IF( pCounter > nParticles ) THEN
            EXIT ParticlePlacement
          END IF
        END DO
      END DO
    END DO
  END DO
END DO ParticlePlacement

! Correct coordinates inside the triclinic simulation box
DO pParticle = 1, nParticles
  CALL MatrixVectorMultiplication( CellLength, pPosition(:,pParticle), pPosition(:,pParticle) )
END DO

! Centralize unit box at origin of coordination system (0,0,0)
pPosition(1,:) = pPosition(1,:) - 0.5D0 * ( BoxLength(1) + BoxLength(4) + BoxLength(7) )
pPosition(2,:) = pPosition(2,:) - 0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) )
pPosition(3,:) = pPosition(3,:) - 0.5D0 * ( BoxLength(3) + BoxLength(6) + BoxLength(9) )

! Define equilibrium positions of all nonconvex particles (centers of mass)
pPositionEq(:,:) = pPosition(:,:)

! Position of the cylindrical petals (centers of mass)
DO pParticle = 1, nParticles
  ! Apply periodic boundary conditions
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,pParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,pParticle) )
  ! Reposition of cylinders
  CALL VectorRotation( cPositionBasis(:,1), pQuaternion(:,pParticle), cRotatedPositionBasis(:,1) )
  cPosition(:,1,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,1) ! First quarter
  CALL VectorRotation( cPositionBasis(:,2), pQuaternion(:,pParticle), cRotatedPositionBasis(:,2) )
  cPosition(:,2,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,2) ! Second quarter
  CALL VectorRotation( cPositionBasis(:,3), pQuaternion(:,pParticle), cRotatedPositionBasis(:,3) )
  cPosition(:,3,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,3) ! Third quarter
  CALL VectorRotation( cPositionBasis(:,4), pQuaternion(:,pParticle), cRotatedPositionBasis(:,4) )
  cPosition(:,4,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,4) ! Fourth quarter
END DO

RETURN

END SUBROUTINE ConfigFloppyBox


! *********************************************************************************************** !
!            This subroutine allocates particles according to the LC configuration                !
! *********************************************************************************************** !
SUBROUTINE ConfigLiquidCrystal(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle  ! Counter (particles)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis ! Relative distance (cylinders)

! Liquid-crystal data
OPEN( Unit= 10, File= "Liquid-crystal/lc.dat", Action= "READ" )
! Skip
READ( 10, * )
READ( 10, * )
READ( 10, * )
! Box length
READ( 10, * ) DummyText, BoxLength
! Box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )
! Skip
READ( 10, * )
! Particle positions
READ( 10, * ) DummyText, pPosition
! Particle quaternions
READ( 10, * ) DummyText, pQuaternion
! Generate orientations along z-direction
DO pParticle = 1, nParticles
  CALL VectorRotation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
END DO
! Generate petals (cylindrical subunits)
DO pParticle = 1, nParticles
  CALL VectorRotation( cPositionBasis(:,1), pQuaternion(:,pParticle), cRotatedPositionBasis(:,1) )
  cPosition(:,1,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,1) ! First quarter
  CALL VectorRotation( cPositionBasis(:,2), pQuaternion(:,pParticle), cRotatedPositionBasis(:,2) )
  cPosition(:,2,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,2) ! Second quarter
  CALL VectorRotation( cPositionBasis(:,3), pQuaternion(:,pParticle), cRotatedPositionBasis(:,3) )
  cPosition(:,3,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,3) ! Third quarter
  CALL VectorRotation( cPositionBasis(:,4), pQuaternion(:,pParticle), cRotatedPositionBasis(:,4) )
  cPosition(:,4,pParticle) = pPosition(:,pParticle) + cRotatedPositionBasis(:,4) ! Fourth quarter
END DO
CLOSE( 10 )

RETURN

END SUBROUTINE ConfigLiquidCrystal

! *********************************************************************************************** !
!          This subroutine creates a file with all particles' positions and orientations          !
! *********************************************************************************************** !
SUBROUTINE OutputConfiguration(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle, cCylinder ! Counters
INTEGER( Kind= Int64 ) :: FrameLEFT            ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT           ! Box frame dimension

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: WritePosition ! Position of particles (center of mass)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Floppy-box structure
IF( ConfigSelectionLogical(1) ) THEN
  OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_initconf_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//"_fb.xyz" )
  WRITE( 10, "(I5)" ) nParticles * 4
  ! Descriptor string
  DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO pParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of fake cylinders
      WritePosition(1) = cPosition(1,cCylinder,pParticle)
      WritePosition(2) = cPosition(2,cCylinder,pParticle)
      WritePosition(3) = cPosition(3,cCylinder,pParticle)
      WRITE( 10, "(G0,10(' ',G0.6))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,pParticle), &
      &                                pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), &
      &                                0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
    END DO
  END DO
  FLUSH( 10 )
  CLOSE( 10 )
ELSE IF( ConfigSelectionLogical(2) ) THEN
  OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_initconf_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//"_lc.xyz" )
  WRITE( 10, "(I5)" ) nParticles * 4
  ! Descriptor string
  DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO pParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of fake cylinders
      WritePosition(1) = cPosition(1,cCylinder,pParticle)
      WritePosition(2) = cPosition(2,cCylinder,pParticle)
      WritePosition(3) = cPosition(3,cCylinder,pParticle)
      WRITE( 10, "(G0,10(' ',G0.6))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,pParticle), &
      &                                pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), &
      &                                0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
    END DO
  END DO
  FLUSH( 10 )
  CLOSE( 10 )
END IF

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 34 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 34 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"INITIAL CONFIGURATION FILE (OVITO)"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(G0)" ) "Initial configuration already available!"
WRITE( *, "(G0)" ) "Check 'Initial_Configuration/OVITO/' folder for more information."
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE OutputConfiguration

END MODULE InitialConfiguration