! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!      This module initialize common variables (number of particles, reduced number density,      !
!       reduced temperature etc.), Monte Carlo parameters (total number of cycles, number of      !
! equilibration cycles etc.), and Block Averaging parameters (maximum/minimum number of blocks).  !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
!                        which results are written out in external files.                         !
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
! Main Reference:                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE VariableInitialization

! Uses two MODULES: global variables and initial configuration
USE GlobalVariables
USE InitialConfiguration

IMPLICIT NONE

! *********************************************************************************************** !
!                                     VARIABLE INITIALIZATION                                     !
!                    Most variables should be first specified in an input file.                   !
!              Some variables are on-the-fly (OTF) and must be specified by the user              !
!                               at the beginning of the simulation.                               !
!                       We provided an example .ini file to guide the user.                       !
!                     Please also check our README file for more information.                     !
! *********************************************************************************************** !
CONTAINS

! *********************************************************************************************** !
!                           INITIALIZATION OF COMMON/GENERAL VARIABLES                            !
! *********************************************************************************************** !
SUBROUTINE GeneralVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: FrameLEFT  ! Box frame dimension (west side)
INTEGER( Kind= Int64 ) :: FrameRIGHT ! Box frame dimension (east side)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: SurfaceArea ! Surface area of the considered molecular geometry

! System properties
OPEN( Unit= 10, File= "ini_system.ini", Action= "READ" )

! Number of particles
READ( 10, * ) DummyText, nParticles
! Reduced number density
READ( 10, * ) DummyText, nDensity
! Reduced spring constant
READ( 10, * ) DummyText, rSpringConstant
! Orientational field stiffness
READ( 10, * ) DummyText, OrientFieldStrength
! Reduced spring constant
READ( 10, * ) DummyText, OrientFieldVector
! Reduced temperature
READ( 10, * ) DummyText, rTemperature
! Calculation type
READ( 10, * ) DummyText, CalculationType

CLOSE( 10 )

! Cross-sectional area
CALL CrossSectionalArea( SurfaceArea )

! Molecular volume
pVolume = SurfaceArea * cLength ! Å³

! Packing fraction
PackingFraction = nDensity * pVolume

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 17 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 17 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"SYSTEM PROPERTIES"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(2G0)" ) "Number of Nonconvex Particles: ", nParticles
WRITE( *, "(2G0)" ) "Number of Cylinders (Petals): ", 4 * nParticles
WRITE( *, "(G0,G0.5,G0)" ) "Cylindrical Diameter (D): ", cDiameter, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Cylindrical Length (L): ", cLength, "Å"
WRITE( *, "(G0,G0.5)" ) "Cylindrical Aspect Ratio (L/D): ", cAspectRatio
WRITE( *, "(G0,G0.5,G0)" ) "Number Density of Particles: ", nDensity, "Å⁻³"
WRITE( *, "(G0,G0.5,G0)" ) "Surface Area (Four-Leaf Clover): ", SurfaceArea, "Å²"
WRITE( *, "(G0,G0.5,G0)" ) "Molecular Volume: ", pVolume, "Å³"
WRITE( *, "(G0,G0.5)" ) "Packing Fraction: ", PackingFraction
WRITE( *, "(G0,G0.5)" ) "Reduced Spring Constant: ", rSpringConstant
WRITE( *, "(G0,G0.5)" ) "Orientational Field Stiffness: ", OrientFieldStrength
WRITE( *, "(G0,2(G0.5,G0),G0.5)" ) "Orientational Field Vector: ", OrientFieldVector(1), ", ", OrientFieldVector(2), ", ", &
&                                                                  OrientFieldVector(3)
WRITE( *, "(G0,G0.5)" ) "Reduced Temperature: ", rTemperature
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 32 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 32 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"INITIAL CONFIGURATION PROPERTIES"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
IF( ConfigName == "FB" ) WRITE( *, "(3G0)" ) "Initial Configuration: Floppy Box (", TRIM( ConfigName ), ")"
IF( ConfigName == "LC" ) WRITE( *, "(3G0)" ) "Initial Configuration: Liquid Crystal (", TRIM( ConfigName ), ")"
WRITE( *, "(2G0)" ) "Arrangement Type: ", cArrangement
WRITE( *, "(G0)" ) " "

! Output file descriptors (Diameter [1])
FileFormat(1) = "(F0.5)"
WRITE( FileDescriptor(1), FileFormat(1) ) cDiameter
! Output file descriptors (cLength [2])
FileFormat(2) = "(F0.5)"
WRITE( FileDescriptor(2), FileFormat(2) ) cLength
! Output file descriptors (Aspect ratio [3])
FileFormat(3) = "(F0.5)"
WRITE( FileDescriptor(3), FileFormat(3) ) cAspectRatio

RETURN

END SUBROUTINE GeneralVariables

! *********************************************************************************************** !
!                            INITIALIZATION OF MONTE CARLO PARAMETERS                             !
! *********************************************************************************************** !
SUBROUTINE MonteCarloVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: FrameLEFT  ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT ! Box frame dimension

! Simulation properties
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )

! Skip line
READ( 10, * ) DummyText, DummyText
! Total number of cycles
READ( 10, * ) DummyText, nSimulationCycles
! Number of equilibration cycles
READ( 10, * ) DummyText, nEquilCycles
! Saving frequency
READ( 10, * ) DummyText, nSavingFrequency
! Adjustment frequency
READ( 10, * ) DummyText, nAdjustmentFrequency
! Maximum translational displacement
READ( 10, * ) DummyText, maxTranslationalDisplc
! Maximum rotational displacement
READ( 10, * ) DummyText, maxRotationalDisplc

CLOSE( 10 )

! Simulation properties
OPEN( Unit= 10, File= "ini_probabilities.ini", Action= "READ" )

! Translational probability
READ( 10, * ) DummyText, TranslationalProbability
! Rotational probability
RotationProbability = 1.D0 - TranslationalProbability

CLOSE( 10 )

! Simulation properties
OPEN( Unit= 10, File= "ini_ratios.ini", Action= "READ" )

! Translational acceptance ratio threshold
READ( 10, * ) DummyText, AcceptanceRatioTranslation
! Rotational acceptance ratio threshold
READ( 10, * ) DummyText, AcceptanceRatioRotation

CLOSE( 10 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 21 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 21 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"SIMULATION PROPERTIES"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(2G0)" ) "Number of Cycles: ", nSimulationCycles
WRITE( *, "(2G0)" ) "Number of Equilibration Cycles: ", nEquilCycles
WRITE( *, "(2G0)" ) "Number of Production Cycles: ", nSimulationCycles - nEquilCycles
WRITE( *, "(3G0)" ) "Saving Frequency: Every ", nSavingFrequency, " Cycle(s)"
WRITE( *, "(3G0)" ) "Adjustment Frequency: Every ", nAdjustmentFrequency, " Cycle(s)"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Translational Displacement: ", maxTranslationalDisplc, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Rotational Displacement: ", maxRotationalDisplc, "rad"
WRITE( *, "(G0,G0.5)" ) "Translational Acceptance Ratio Threshold: ", AcceptanceRatioTranslation
WRITE( *, "(G0,G0.5)" ) "Rotational Acceptance Ratio Threshold: ", AcceptanceRatioRotation
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Translational Displacement: ", TranslationalProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Rotational Displacement: ", RotationProbability * 100.D0, "%"
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE MonteCarloVariables

! *********************************************************************************************** !
!                           INITIALIZATION OF INQUIRY/CONTROL VARIABLES                           !
! *********************************************************************************************** !
SUBROUTINE ControlVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: FrameLEFT  ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT ! Box frame dimension

! Control variables
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )

! Skip
READ( 10 , * ) DummyText, DummyText
READ( 10 , * ) DummyText, DummyNumber
READ( 10 , * ) DummyText, DummyNumber
READ( 10 , * ) DummyText, DummyNumber
READ( 10 , * ) DummyText, DummyNumber
READ( 10 , * ) DummyText, DummyNumber
READ( 10 , * ) DummyText, DummyNumber

! Trajectory inquiry
READ( 10 , * ) DummyText, TrajectoryInquiry
CALL ToUpper( TrajectoryInquiry, LEN_TRIM( TrajectoryInquiry ), TrajectoryInquiry )
IF( TrajectoryInquiry == "Y" ) THEN
  TrajectoryLogical = .TRUE.
ELSE IF( TrajectoryInquiry == "N" ) THEN
  TrajectoryLogical = .FALSE.
END IF

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 18 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 18 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"CONTROL PROPERTIES"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(3G0)" ) "Print Trajectory? [", TrajectoryInquiry, "]"
IF( CalculationType == 0 ) THEN
  WRITE( *, "(2G0)" ) "Calculation Type (Frenkel-Ladd): A", CalculationType
ELSE IF( CalculationType == 1 ) THEN
  WRITE( *, "(2G0)" ) "Calculation Type (Frenkel-Ladd): ΔA", CalculationType
ELSE IF( CalculationType == 2 ) THEN
  WRITE( *, "(2G0)" ) "Calculation Type (Frenkel-Ladd): ΔA", CalculationType
ELSE IF( CalculationType == 3 ) THEN
  WRITE( *, "(G0)" ) "Calculation Type (Frenkel-Mulder): ΔAfield"
END IF
WRITE( *, "(G0)" ) "Potential Type: Hard-Core Purely Repulsive Potential"
WRITE( *, "(G0)" ) " "

CLOSE( 10 )

RETURN

END SUBROUTINE ControlVariables

END MODULE VariableInitialization