! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!  This program determines the free energy of a solid/nematic liquid crystal through free-energy  !
!  methods, including the Frenkel-Ladd method (C. Vega, E. Sanz, Abascal J. & Noya E., 2008) and  !
!                   the Frenkel-Mulder method (D. Frenkel, B. M. Mulder, 1985).                   !
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
! Main References:                        C. Vega, S. Lago                                        !
!                                 Computers Chem. 18, 55-59 (1993)                                !
!                                DOI: 10.1016/0097-8485(94)80023-5                                !
!                             --------------------------------------                              !
!                            A. G. Orellana, E. Romani, C. de Michele                             !
!                                  Eur. Phys. J. E 41, 10 (2018)                                  !
!                                 DOI: 10.1140/epje/i2018-11657-0                                 !
!                             --------------------------------------                              !
!                J. T. Lopes, F. Romano, E. Grelet, L. F. M. Franco, A. Giacometti                !
!                                J. Chem. Phys. 154, 104902 (2021)                                !
!                                     DOI: 10.1063/5.0040942                                      !
!                             --------------------------------------                              !
!                                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                         C. Vega, E. Sanz, J. L. F. Abascal, E. G. Noya                          !
!                           J. Phys.: Condens. Matter 20, 153101 (2008)                           !
!                               DOI: 10.1088/0953-8984/20/15/153101                               !
!                             --------------------------------------                              !
!                                    D. Frenkel, B. M. Mulder                                     !
!                               J. Mol. Phys. 55, 1171-1192 (1985)                                !
!                                 DOI: 10.1080/00268978500101971                                  !
!                             --------------------------------------                              !
!                                      N. Metropolis et al.                                       !
!                                 J. Chem. Phys. 21, 1087 (1953)                                  !
!                                     DOI: 10.1063/1.1699114                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
PROGRAM FreeEnergyMethods

! Uses five modules: global variables, variable initialization, initial configuration, directory creator, and thermodynamic path
USE GlobalVariables
USE VariableInitialization
USE InitialConfiguration
USE Folders
USE ThermodynamicStep

IMPLICIT NONE

! Get the maximum number of threads
#ifdef _OPENMP
  nThreads = OMP_GET_MAX_THREADS(  )
#else
  nThreads = 1
#endif

! *********************************************************************************************** !
! Molecular configuration selection (see 'InitialConfiguration' module)                           !
! *********************************************************************************************** !
CALL ConfigurationSelection(  )

! *********************************************************************************************** !
! Cross-section arrangement selection (see 'InitialConfiguration' module)                         !
! *********************************************************************************************** !
CALL GeometrySelection(  )

! Unrotated reference position (cylinders)
IF( GeometrySelectionLogical(1) ) THEN ! GEOMETRY [1]
  ! First quarter
  cPositionBasis(1,1) = 0.25D0 * cDiameter
  cPositionBasis(2,1) = 0.25D0 * cDiameter
  cPositionBasis(3,1) = 0.D0
  ! Second quarter
  cPositionBasis(1,2) = -0.25D0 * cDiameter
  cPositionBasis(2,2) = 0.25D0 * cDiameter
  cPositionBasis(3,2) = 0.D0
  ! Third quarter
  cPositionBasis(1,3) = -0.25D0 * cDiameter
  cPositionBasis(2,3) = -0.25D0 * cDiameter
  cPositionBasis(3,3) = 0.D0
  ! Fourth quarter
  cPositionBasis(1,4) = 0.25D0 * cDiameter
  cPositionBasis(2,4) = -0.25D0 * cDiameter
  cPositionBasis(3,4) = 0.D0
ELSE IF( GeometrySelectionLogical(2) ) THEN ! GEOMETRY [2]
  ! First quarter
  cPositionBasis(1,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,1) = 0.D0
  ! Second quarter
  cPositionBasis(1,2) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,2) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,2) = 0.D0
  ! Third quarter
  cPositionBasis(1,3) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,3) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,3) = 0.D0
  ! Fourth quarter
  cPositionBasis(1,4) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,4) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,4) = 0.D0
END IF

! *********************************************************************************************** !
! CPU Clock                                                                                       !
! *********************************************************************************************** !
CALL Date_and_Time( Values= DateTime )

! *********************************************************************************************** !
! Initialization of general variables (see 'VariableInitialization' module)                       !
! *********************************************************************************************** !
CALL GeneralVariables(  )

! Allocation
ALLOCATE( pQuaternion(0:3,nParticles) )
ALLOCATE( pPosition(3,nParticles) )
ALLOCATE( pOrientation(3,nParticles) )
ALLOCATE( cPosition(3,4,nParticles) )
ALLOCATE( pQuaternionMC(0:3,nParticles) )
ALLOCATE( pPositionMC(3,nParticles) )
ALLOCATE( pOrientationMC(3,nParticles) )
ALLOCATE( cPositionMC(3,4,nParticles) )

! *********************************************************************************************** !
! Initialization of Monte Carlo parameters (see 'VariableInitialization' module)                  !
! *********************************************************************************************** !
CALL MonteCarloVariables(  )

! *********************************************************************************************** !
! Initialization of Inquiry/Control variables (see 'VariableInitialization' module)               !
! *********************************************************************************************** !
CALL ControlVariables(  )

! *********************************************************************************************** !
! Initial configuration folder (see 'Folders' module)                                             !
! *********************************************************************************************** !
CALL InitialConfigFolder(  )

! *********************************************************************************************** !
! Create output directories (see 'Folders' module)                                                !
! *********************************************************************************************** !
CALL MainFolders(  )

! *********************************************************************************************** !
! Create date subfolders (see 'Folders' module)                                                   !
! *********************************************************************************************** !
CALL DateFolders(  )

! *********************************************************************************************** !
! Initial configuration (see 'InitialConfiguration' module)                                       !
! *********************************************************************************************** !
IF( ConfigSelectionLogical(1) ) THEN ! Calls 'ConfigFloppyBox' subroutine if the user is computing the free energy of a solid
  CALL ConfigFloppyBox(  )
ELSE IF( ConfigSelectionLogical(2) ) THEN ! Calls 'ConfigLiquidCrystal' subroutine if the user is computing the free energy of a liquid crystal
  CALL ConfigLiquidCrystal(  )
END IF

! Summary
WRITE( *, "(2G0)" ) "Number of threads: ", nThreads
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Proceed? [Y/N]"
READ( *, * ) DummyText
WRITE( *, "(G0)" ) " "
CALL ToUpper( DummyText, LEN_TRIM( DummyText ), DummyText )
IF( DummyText == "N" ) THEN
  STOP
END IF

! *********************************************************************************************** !
! Initial configuration file (see 'InitialConfiguration' module)                                  !
! *********************************************************************************************** !
CALL OutputConfiguration(  )

! Cutoff diameter (non-convex body)
IF( GeometrySelectionLogical(1) ) THEN
  CutoffSphere = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter + cLength )
  SquaredCutoffSphere = CutoffSphere * CutoffSphere
  CutoffSpherocylinder = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter )
  SquaredCutoffSpherocylinder = CutoffSpherocylinder * CutoffSpherocylinder
ELSE IF( GeometrySelectionLogical(2) ) THEN
  CutoffSphere = ( 2.D0 * cDiameter + cLength )
  SquaredCutoffSphere = CutoffSphere * CutoffSphere
  CutoffSpherocylinder = 2.D0 * cDiameter
  SquaredCutoffSpherocylinder = CutoffSpherocylinder * CutoffSpherocylinder
END IF

! Cutoff diameter (cylinders)
cCutoffSphere = ( cDiameter + cLength )
cSquaredCutoffSphere = cCutoffSphere * cCutoffSphere

! *********************************************************************************************** !
! Objective                                                                                       !
! *********************************************************************************************** !
IF( CalculationType == 0 ) THEN
  ! Calculate A0 (free energy change between an unconstrained solid and the solid with fixed CM + ideal contributions)
  CALL Compute_A0(  )
ELSE IF( CalculationType == 1 ) THEN
  ! Calculate ΔA1 (free energy change between an interacting EC and a non-interacting EC)
  CALL Compute_A1(  )
ELSE IF( CalculationType == 2 ) THEN
  ! Calculate ΔA2 (free energy change between the solid and the interacting EC)
  CALL Compute_A2(  )
ELSE IF( CalculationType == 3 ) THEN
  ! Calculate ΔAfield (free energy change between the externally oriented fluid and the liquid crystal)
  CALL Compute_AField(  )
END IF

! Status
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Finished!"

END PROGRAM FreeEnergyMethods