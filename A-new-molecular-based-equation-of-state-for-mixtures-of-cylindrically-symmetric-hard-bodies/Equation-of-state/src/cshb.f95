! ************************************************************************************************ !
!                          CSHB Equation of State (Square-Well Potential)                          !
! ************************************************************************************************ !
!        This program is used to determine thermodynamic properties and phase equilibria of        !
!                         substances through the CSHB-SW equation of state                         !
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
PROGRAM CSHB_EOS

USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************* !
! INTEGER VARIABLES                                                                     !
! ************************************************************************************* !
INTEGER( Kind= Int64 ) :: iComponent, jComponent, cComponent ! Counter (components)
INTEGER( Kind= Int64 ) :: Path                               ! Path selection

! ************************************************************************************* !
! REAL VARIABLES                                                                        !
! ************************************************************************************* !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: mFraction ! Molar fraction

! ************************************************************************************* !
! CHARACTER STRINGS                                                                     !
! ************************************************************************************* !
CHARACTER( Len= 01 ) :: Dummy ! Dummy

! Get date and time
CALL Get_Date_and_Time(  )

! Read the specification file
CALL Read_Global_Specification(  )

! Allocate shared arrays
CALL Allocation_Arrays(  )

! Allocate private arrays
ALLOCATE( mFraction(nComponents) )

! Component selection
DO cComponent = 1, nComponents
  WRITE( *, "(3G0)" ) "Enter name of component #", cComponent, ": "
  READ( *, * ) cMoleculeName(cComponent)
  CALL ToUpper( cMoleculeName(cComponent), LEN_TRIM( cMoleculeName(cComponent) ), cMoleculeName(cComponent) )
  IF( nComponents > 1 .AND. cComponent > 1 ) THEN
    DO iComponent = 1, nComponents - 1
      DO jComponent = iComponent + 1, nComponents
        IF( TRIM( cMoleculeName(iComponent) ) == TRIM( cMoleculeName(jComponent) ) ) THEN
          WRITE( *, "(G0)" ) "Components are the same. Exiting..."
          CALL Exit(  )
        END IF
      END DO
    END DO
  END IF
END DO

! Geometry selection
GeometrySpecification = .FALSE.
DO cComponent = 1, nComponents
  IF( cGeometry(cComponent) == 1 ) THEN
    GeometrySpecification(cComponent,1) = .TRUE.
  ELSE IF( cGeometry(cComponent) == 2 ) THEN
    GeometrySpecification(cComponent,2) = .TRUE.
  ELSE IF( cGeometry(cComponent) == 3 ) THEN
    GeometrySpecification(cComponent,3) = .TRUE.
  END IF
END DO

! Component properties
DO cComponent = 1, nComponents
  CALL Properties( cComponent, cMoleculeName(cComponent) )
  CALL Specific_Heat_Parameters( cComponent, cMoleculeName(cComponent) )
  CALL Structural_Formula( cComponent, cMoleculeName(cComponent) )
END DO

! Mixture properties
CALL Read_Mixture_Specification( mFraction )

! Mixing and combining rules
CALL Combining_Rules(  )
CALL Combined_Properties(  )
CALL Isihara_Hadwiger_Theorem(  )

! Pure component analysis
IF( nComponents == 1 ) THEN
  ! Status
  WRITE( *, "(G0)" ) "Pure Component Analysis"
  WRITE( *, "(G0)" ) " "
  ! Read the path file
  OPEN( Unit= 10, File= "path.ini", Action= "READ" )
  READ( 10, * ) Dummy, Path
  IF( Path < 1 .OR. Path > 3 ) THEN
    WRITE( *, "(2G0)" ) "Error: the thermodynamic path must be 1 (thermodynamic properties), 2 (coexistence curve), or 3 ", &
    &                   "(optimization). Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 10 )
  ! Paths
  IF( Path == 1 ) THEN
    CALL Thermodynamic_Properties( .TRUE. )
  ELSE IF( Path == 2 ) THEN
    CALL Coexistence_Curve(  )
  ELSE IF( Path == 3 ) THEN
    CALL Simplex_Optimization( .TRUE. )
  END IF
! Multicomponent analysis
ELSE
  ! Status
  WRITE( *, "(G0)" ) "Multicomponent Analysis"
  WRITE( *, "(G0)" ) " "
  ! Read the path file
  OPEN( Unit= 10, File= "path.ini", Action= "READ" )
  READ( 10, * ) Dummy, Path
  IF( Path < 1 .OR. Path > 3 ) THEN
    WRITE( *, "(2G0)" ) "Error: the thermodynamic path must be 1 (thermodynamic properties), 2 (coexistence curve), or 3 ", &
    &                   "(optimization). Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 10 )
  ! Paths
  IF( Path == 1 ) THEN
    CALL Thermodynamic_Properties( .FALSE. )
  ELSE IF( Path == 2 ) THEN
    CALL Binary_Vapor_Liquid_Equilibrium(  )
  ELSE IF( Path == 3 ) THEN
    CALL Simplex_Optimization( .FALSE. )
  END IF
END IF

END PROGRAM CSHB_EOS