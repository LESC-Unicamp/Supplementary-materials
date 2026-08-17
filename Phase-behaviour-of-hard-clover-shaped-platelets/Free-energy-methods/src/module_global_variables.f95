! ############################################################################################### !
!                                     EINSTEIN CRYSTAL METHOD                                     !
!           This module defines the variables used by the main program and most of the            !
!         subroutines and functions. A brief description is presented for each variable.          !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                        April 11th, 2023                                         !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE GlobalVariables

! Uses two intrinsic modules: Int64, Real64, and Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Int64, Real64, Output_Unit

! OpenMP API
#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nParticles      ! Number of particles
INTEGER( Kind= Int64 ) :: nThreads        ! Number of threads (OpenMP)
INTEGER( Kind= Int64 ) :: cArrangement    ! Cross-sectional geometry
INTEGER( Kind= Int64 ) :: CalculationType ! Calculation type (0 = A0, 1 = ΔA1, 2 = ΔA2, and 3 = ΔAfield)

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nSimulationCycles    ! Total number of cycles
INTEGER( Kind= Int64 ) :: nEquilCycles         ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nSavingFrequency     ! Saving frequency
INTEGER( Kind= Int64 ) :: nAdjustmentFrequency ! Adjustment frequency

! *********************************************************************************************** !
! REAL VARIABLES (GENERAL)                                                                        !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: BoxVolume                   ! Volume of simulation box
REAL( Kind= Real64 )                    :: nDensity                    ! Number density
REAL( Kind= Real64 )                    :: PackingFraction             ! Packing fraction
REAL( Kind= Real64 )                    :: rTemperature                ! Reduced temperature
REAL( Kind= Real64 )                    :: pVolume                     ! Molecular volume
REAL( Kind= Real64 )                    :: cAspectRatio                ! Length-to-diameter (L/D) aspect ratio
REAL( Kind= Real64 )                    :: cDiameter                   ! Diameter (cylinder)
REAL( Kind= Real64 )                    :: cLength                     ! Length (cylinder)
REAL( Kind= Real64 )                    :: rSpringConstant             ! Reduced spring constant
REAL( Kind= Real64 )                    :: RandomNumber                ! Pseudorandom number
REAL( Kind= Real64 )                    :: TranslationalProbability    ! Probability of movement (translation)
REAL( Kind= Real64 )                    :: RotationProbability         ! Probability of movement (rotation)
REAL( Kind= Real64 )                    :: DummyNumber                 ! Dummy (number)
REAL( Kind= Real64 )                    :: OrientFieldStrength         ! Strength of the orientational field
REAL( Kind= Real64 )                    :: CutoffSphere                ! Cutoff diameter of the circumscribing sphere (particle)
REAL( Kind= Real64 )                    :: cCutoffSphere               ! Cutoff diameter of the circumscribing sphere (cylinder)
REAL( Kind= Real64 )                    :: SquaredCutoffSphere         ! Squared cutoff diameter of the circumscribing sphere (particle)
REAL( Kind= Real64 )                    :: cSquaredCutoffSphere        ! Squared cutoff diameter of the circumscribing sphere (cylinder)
REAL( Kind= Real64 )                    :: CutoffSpherocylinder        ! Cutoff diameter of the circumscribing spherocylinder (particle)
REAL( Kind= Real64 )                    :: SquaredCutoffSpherocylinder ! Squared cutoff diameter of the circumscribing spherocylinder (particle)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: OrientFieldVector           ! Orientation vector of the orientational field
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLength                   ! Length (x,y,z) of simulation box (triclinic)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverse            ! Inverse of length (x,y,z) of simulation box (triclinic)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bPositionEq                 ! Equilibrium position of the center of mass of the system
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cPositionBasis              ! Body-fixed position

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS)                                                         !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: maxTranslationalDisplc     ! User maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: maxRotationalDisplc        ! User maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: AcceptanceRatioTranslation ! Acceptance ratio threshold (Translation)
REAL( Kind= Real64 ) :: AcceptanceRatioRotation    ! Acceptance ratio threshold (Rotation)

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS, ALLOCATABLE)                                            !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPosition, pPositionMC       ! Position array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pOrientation, pOrientationMC ! Orientation array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternion, pQuaternionMC   ! Quaternion array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPositionCell                ! Equilibrium position of particles in the unit cell
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPositionEq                  ! Equilibrium position of particles in the crystal
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternionCell              ! Equilibrium orientation (quaternion) of particles in the unit cell
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternionEq                ! Equilibrium orientation (quaternion) of particles in the crystal
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: cPosition, cPositionMC       ! Position array (cylinders)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: pOrientationCell             ! Equilibrium orientation of particles in the unit cell [three orthogonal vectors]
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: pOrientationEq               ! Equilibrium orientation of particles in the crystal [three orthogonal vectors]

! *********************************************************************************************** !
! REAL PARAMETERS (CONSTANTS)                                                                     !
! *********************************************************************************************** !
REAL( Kind= Real64 ), PARAMETER                 :: cPi = 4.D0 * DATAN( 1.D0 ) ! π
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: xAxis = [1.D0, 0.D0, 0.D0] ! Body-fixed axis of rotation (x)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: yAxis = [0.D0, 1.D0, 0.D0] ! Body-fixed axis of rotation (y)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: zAxis = [0.D0, 0.D0, 1.D0] ! Body-fixed axis of rotation (z)

! *********************************************************************************************** !
! CHARACTER STRINGS (GENERAL)                                                                     !
! *********************************************************************************************** !
CHARACTER( LEN= 01 ) :: TrajectoryInquiry ! Trajectory output inquiry
CHARACTER( LEN= 03 ) :: ConfigName        ! Molecular configuration inquiry
CHARACTER( LEN= 01 ) :: DummyText         ! Dummy (character)

! *********************************************************************************************** !
! CHARACTER STRINGS (FILE/FOLDER ORGANIZER)                                                       !
! *********************************************************************************************** !
CHARACTER( LEN= 10 ) :: HourDescriptor                 ! Descriptor for output file (date and time)
CHARACTER( LEN= 10 ) :: DateDescriptor                 ! Descriptor for output folder (date and time)
CHARACTER( LEN= 32 ) :: HourFormat                     ! String format for output file (date and time)
CHARACTER( LEN= 32 ) :: DateFormat                     ! String format for output folder (date and time)
CHARACTER( LEN= 32 ), DIMENSION( 4 ) :: FileFormat     ! String format for output file
CHARACTER( LEN= 32 ), DIMENSION( 4 ) :: FileDescriptor ! Descriptor for output file

! *********************************************************************************************** !
! CHARACTER PARAMETERS (BOX FRAMES)                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 3 ), PARAMETER :: CH_HS = "═" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VS = "║" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_UL = "╔" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_BL = "╚" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_UR = "╗" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_BR = "╝" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VL = "╠" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VR = "╣" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SH_VL = "╟" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SH_VR = "╢" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_HS = "─" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VS = "│" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VR = "├" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VL = "┤" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_UL = "┌" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_BL = "└" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_UR = "┐" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_BR = "┘" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: C_FUL = "█" ! Box drawing symbol

! *********************************************************************************************** !
! LOGICAL VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
LOGICAL                 :: TrajectoryLogical        ! Trajectory output selection
LOGICAL, DIMENSION( 2 ) :: ConfigSelectionLogical   ! Molecular configuration selection
LOGICAL, DIMENSION( 2 ) :: GeometrySelectionLogical ! Cross-sectional geometry selection

! *********************************************************************************************** !
! LOGICAL VARIABLES (FILE/FOLDER ORGANIZER)                                                       !
! *********************************************************************************************** !
LOGICAL                 :: FileExistLogical       ! Checks whether a file exists or not
LOGICAL, DIMENSION( 5 ) :: FolderExistLogical     ! Checks whether folder exists or not
LOGICAL, DIMENSION( 4 ) :: DateFolderExistLogical ! Checks whether date folders exist or not
LOGICAL, DIMENSION( 3 ) :: SubfolderExistLogical  ! Checks whether subfolder exists or not

END MODULE GlobalVariables