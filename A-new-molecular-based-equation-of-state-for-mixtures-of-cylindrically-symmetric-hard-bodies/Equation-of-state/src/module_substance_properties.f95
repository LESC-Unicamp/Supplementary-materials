! ************************************************************************************************ !
!                                 DATABASE OF CHEMICAL SUBSTANCES                                  !
! ************************************************************************************************ !
!       This modules is a database for the thermodynamic properties of all molecules used in       !
!                         the main program and its underlying subroutines.                         !
! ************************************************************************************************ !
! => AUTHOR:     Nathan Barros de Souza                                                            !
! => E-MAIL:     n264179@dac.unicamp.br                                                            !
! => SUPERVISOR: Luís Fernando Mercier Franco                                                      !
! ************************************************************************************************ !
! Main References:         B. E. Poling, J. M. Prausnitz, J. P. O'Connell                          !
!           "The Properties of Gases and Liquids", 5th Ed., McGraw-Hill Education (2001)           !
!                                                                                                  !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !
MODULE Substances

IMPLICIT NONE

CONTAINS

! ************************************************************************************************ !
!                         Defines all molecular properties within database                         !
! ************************************************************************************************ !
SUBROUTINE Properties( cComponent, cName )

! Use some variables of the global variables module
USE GlobalVariables, ONLY: Int64, Real64, ZhangCorrectionLogical, HigherOrderTPTLogical, GeometrySpecification, cDiameter, &
&                          cWellDepth, cPotentialRange, cAspectRatio, cLength, cMolarMass, cBoltzmann, cDiameterSphere, &
&                          PotentialTypeLogical, fDiameter, fLength, fAspectRatio, PYHCBCorrectionLogical, &
&                          ReferenceBoublikLogical, UseA1ForA2Logical

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index
INTEGER( Kind= Int64 ) :: DataSet    ! Data set index

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( LEN= * ) :: cName ! Component name

! Database of substances
IF( .NOT. HigherOrderTPTLogical ) THEN ! No higher-order TPT
  IF( ZhangCorrectionLogical ) THEN ! Zhang correction
    IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-Revolution
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 4.403281697D-10 ! Meter
          cWellDepth(cComponent) = 147.2947075D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.548285565D0
          cAspectRatio(cComponent) = 0.5081012779D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 4.272190676D-10 ! Meter
            cWellDepth(cComponent) = 150.2692255D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.535518934D0
            cAspectRatio(cComponent) = 0.5544676698D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 4.183481033D-10 ! Meter
                cWellDepth(cComponent) = 151.7664481D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.4965811130D0
                cAspectRatio(cComponent) = 0.6193525806D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 4.056247308D-10 ! Meter
                cWellDepth(cComponent) = 150.4977130D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.5019368642D0
                cAspectRatio(cComponent) = 0.6784056846D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.747134326D-10 ! Meter
          cWellDepth(cComponent) = 288.4317847D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.495754227D0
          cAspectRatio(cComponent) = 3.003571151D0
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.945059597D-10 ! Meter
                cWellDepth(cComponent) = 335.0299654D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3534845706D0
                cAspectRatio(cComponent) = 2.572880165D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.906D-10 ! Meter
                cWellDepth(cComponent) = 321.88D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3851D0
                cAspectRatio(cComponent) = 2.5263D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.863377340D-10 ! Meter
          cWellDepth(cComponent) = 394.6070878D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.459307611D0
          cAspectRatio(cComponent) = 3.583416681D0
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 6.130481958D-10 ! Meter
          cWellDepth(cComponent) = 350.1836234D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.590582597D0
          cAspectRatio(cComponent) = 0.2707255997D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 5.774975148D-10 ! Meter
            cWellDepth(cComponent) = 386.4550657D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.526703881D0
            cAspectRatio(cComponent) = 0.3216672936D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 6.168657359D-10 ! Meter
                cWellDepth(cComponent) = 360.9796916D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.5288141247D0
                cAspectRatio(cComponent) = 0.2216660137D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 3.023473498D-10 ! Meter
          cWellDepth(cComponent) = 478.9976412D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.448833512D0
          cAspectRatio(cComponent) = 3.843216605D0
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "CC4H6" .OR. cName == "CYCLOBUTENE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 6.472774661D-10 ! Meter
          cWellDepth(cComponent) = 483.9411301D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.433521703D0
          cAspectRatio(cComponent) = 0.3180120864D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 6.118970316D-10 ! Meter
            cWellDepth(cComponent) = 497.6802988D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.420091448D0
            cAspectRatio(cComponent) = 0.3688516981D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 5.901805471D-10 ! Meter
                cWellDepth(cComponent) = 560.0162756D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3181498001D0
                cAspectRatio(cComponent) = 0.4362376096D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 54.0904D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.077493711D-10 ! Meter
          cWellDepth(cComponent) = 568.6778470D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.438551107D0
          cAspectRatio(cComponent) = 4.357452112D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "CC5H10" .OR. cName == "CYCLOPENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 7.422439757D-10 ! Meter
          cWellDepth(cComponent) = 579.9444699D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.437968982D0
          cAspectRatio(cComponent) = 0.2690782374D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 6.932003976D-10 ! Meter
            cWellDepth(cComponent) = 604.8122650D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.414142451D0
            cAspectRatio(cComponent) = 0.3237819866D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 6.711667376D-10 ! Meter
                cWellDepth(cComponent) = 714.5630385D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2687082562D0
                cAspectRatio(cComponent) = 0.4213309191D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 6.580226036D-10 ! Meter
                cWellDepth(cComponent) = 710.5776168D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2825046921D0
                cAspectRatio(cComponent) = 0.4194285511D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 70.1329D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.161619076D-10 ! Meter
          cWellDepth(cComponent) = 652.8148108D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.426627769D0
          cAspectRatio(cComponent) = 4.715829648D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "CC6H12" .OR. cName == "CYCLOHEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 7.946338336D-10 ! Meter
          cWellDepth(cComponent) = 611.6512416D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.469103854D0
          cAspectRatio(cComponent) = 0.2519825148D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 7.408283682D-10 ! Meter
            cWellDepth(cComponent) = 657.1000076D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.425455934D0
            cAspectRatio(cComponent) = 0.3075518189D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 7.147195059D-10 ! Meter
                cWellDepth(cComponent) = 800.0371396D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2596775608D0
                cAspectRatio(cComponent) = 0.4140711888D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 7.015228159D-10 ! Meter
                cWellDepth(cComponent) = 805.1672700D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2698931868D0
                cAspectRatio(cComponent) = 0.4094789859D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 84.1595D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.234836636D-10 ! Meter
          cWellDepth(cComponent) = 745.8376783D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.407294392D0
          cAspectRatio(cComponent) = 5.079017782D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.280890240D-10 ! Meter
          cWellDepth(cComponent) = 817.4228653D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.406290910D0
          cAspectRatio(cComponent) = 5.459917987D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.609901851D-10 ! Meter
          cWellDepth(cComponent) = 242.1442879D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456620486D0
          cAspectRatio(cComponent) = 3.428447416D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 5.432068137D-10 ! Meter
                cWellDepth(cComponent) = 300.2313024D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2846324257D0
                cAspectRatio(cComponent) = 0.4348294084D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 5.322320835D-10 ! Meter
                cWellDepth(cComponent) = 296.9991638D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2982197027D0
                cAspectRatio(cComponent) = 0.4401465540D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.782917475D-10 ! Meter
          cWellDepth(cComponent) = 339.9787441D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456218486D0
          cAspectRatio(cComponent) = 4.221148503D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.265211312D-10 ! Meter
                cWellDepth(cComponent) = 301.0997057D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.4199204427D0
                cAspectRatio(cComponent) = 5.917956896D0
              ELSE
                cDiameter(cComponent) = 4.196989173D-10 ! Meter
                cWellDepth(cComponent) = 149.0056363D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.5060090704D0
                cAspectRatio(cComponent) = 0.6096028891D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.034941503D-10 ! Meter
          cWellDepth(cComponent) = 467.6038538D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.387546224D0
          cAspectRatio(cComponent) = 4.547682484D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN ! 15K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 1.981688263D-10 ! Meter
          cWellDepth(cComponent) = 313.3355273D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.543113211D0
          cAspectRatio(cComponent) = 4.529111221D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.105838894D-10 ! Meter
            cWellDepth(cComponent) = 355.1313528D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.468358757D0
            cAspectRatio(cComponent) = 3.722125441D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 1.891663176D-10 ! Meter
                cWellDepth(cComponent) = 376.2028497D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3112491604D0
                cAspectRatio(cComponent) = 5.118009521D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 1.802755694D-10 ! Meter
                cWellDepth(cComponent) = 382.8694021D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3414604231D0
                cAspectRatio(cComponent) = 5.076434636D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(2G0)" ) "There are two set of parameters for BENZENE. Set [1] has an aspect ratio of a (DISK) and ", &
          &                   "set [2] has an aspect ratio of a (ROD). Select the dataset for BENZENE..."
          READ( *, * ) DataSet
          IF( DataSet == 1 ) THEN
            cDiameter(cComponent) = 7.488301839D-10 ! Meter
            cWellDepth(cComponent) = 617.0974480D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.478269644D0
            cAspectRatio(cComponent) = 0.2472845373D0
            IF( ReferenceBoublikLogical ) THEN
              cDiameter(cComponent) = 6.970775761D-10 ! Meter
              cWellDepth(cComponent) = 663.5793896D0 * cBoltzmann ! Joule
              cPotentialRange(cComponent) = 1.433196057D0
              cAspectRatio(cComponent) = 0.3031214390D0
            END IF
          ELSE IF( DataSet == 2 ) THEN
            cDiameter(cComponent) = 2.950166266D-10 ! Meter
            cWellDepth(cComponent) = 617.0961486D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.478272399D0
            cAspectRatio(cComponent) = 4.043954647D0
          ELSE
            WRITE( *, "(G0)" ) "Invalid dataset. Exiting..."
            CALL EXIT(  )
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 6.729145283D-10 ! Meter
                cWellDepth(cComponent) = 812.2918957D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2622046799D0
                cAspectRatio(cComponent) = 0.4080595317D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 6.606784376D-10 ! Meter
                cWellDepth(cComponent) = 818.6350302D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2728012965D0
                cAspectRatio(cComponent) = 0.4015758656D0
              ELSE
                cDiameter(cComponent) = 6.742501360D-10 ! Meter
                cWellDepth(cComponent) = 806.2726915D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2600107082D0
                cAspectRatio(cComponent) = 0.4014392670D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(2G0)" ) "There are two set of parameters for TOLUENE. Set [1] has an aspect ratio of 0.247 (DISK) and ", &
          &                   "set [2] has an aspect ratio of 4.044 (ROD). Select the dataset for TOLUENE..."
          READ( *, * ) DataSet
          IF( DataSet == 1 ) THEN
            cDiameter(cComponent) = 8.355733571D-10 ! Meter
            cWellDepth(cComponent) = 730.0737250D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441922903D0
            cAspectRatio(cComponent) = 0.2167680151D0
          ELSE IF( DataSet == 2 ) THEN
            cDiameter(cComponent) = 3.015203752D-10 ! Meter
            cWellDepth(cComponent) = 730.0760864D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441919167D0
            cAspectRatio(cComponent) = 4.613187985D0
          ELSE
            WRITE( *, "(G0)" ) "Invalid dataset. Exiting..."
            CALL EXIT(  )
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.715556100D-10 ! Meter
          cWellDepth(cComponent) = 265.6128293D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.485946220D0
          cAspectRatio(cComponent) = 2.818480073D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.815565894D-10 ! Meter
            cWellDepth(cComponent) = 273.1640212D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.471849482D0
            cAspectRatio(cComponent) = 2.485634591D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.923135083D-10 ! Meter
                cWellDepth(cComponent) = 301.2509436D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3613358245D0
                cAspectRatio(cComponent) = 2.386908296D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.864525366D-10 ! Meter
                cWellDepth(cComponent) = 292.1490733D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3891050584D0
                cAspectRatio(cComponent) = 2.395955435D0
              ELSE
                cDiameter(cComponent) = 2.471353664D-10 ! Meter
                cWellDepth(cComponent) = 286.2777364D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3841311542D0
                cAspectRatio(cComponent) = 3.447277352D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature and 32K above the triple point
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.797042763D-10 ! Meter
          cWellDepth(cComponent) = 374.9795147D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.477596759D0
          cAspectRatio(cComponent) = 3.489080483D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.940736990D-10 ! Meter
            cWellDepth(cComponent) = 390.4138563D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.454300892D0
            cAspectRatio(cComponent) = 2.946129829D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 3.205276730D-10 ! Meter
                cWellDepth(cComponent) = 457.3562288D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3036747595D0
                cAspectRatio(cComponent) = 2.569527162D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.112028068D-10 ! Meter
                cWellDepth(cComponent) = 449.9449964D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3249565074D0
                cAspectRatio(cComponent) = 2.606840019D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN ! 0K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.070753305D-10 ! Meter
          cWellDepth(cComponent) = 18.14468259D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.689038120D0
          cAspectRatio(cComponent) = 0.9996606604D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.645846473D-10 ! Meter
          cWellDepth(cComponent) = 147.2942652D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.548287342D0
          cAspectRatio(cComponent) = 0.8946617139D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.740367472D-10 ! Meter
            cWellDepth(cComponent) = 150.2692650D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.535518637D0
            cAspectRatio(cComponent) = 0.7339218911D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 3.607106480D-10 ! Meter
                cWellDepth(cComponent) = 155.0677133D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.4777699504D0
                cAspectRatio(cComponent) = 0.9290674701D-3
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.601262644D-10 ! Meter
                cWellDepth(cComponent) = 152.9261372D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.4843734545D0
                cAspectRatio(cComponent) = 0.001242421802D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.511336336D-10 ! Meter
          cWellDepth(cComponent) = 288.4319412D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.495753650D0
          cAspectRatio(cComponent) = 1.954363177D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.737190393D-10 ! Meter
                cWellDepth(cComponent) = 320.3800633D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3520672550D0
                cAspectRatio(cComponent) = 1.488767920D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.719D-10 ! Meter
                cWellDepth(cComponent) = 309.55D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3816D0
                cAspectRatio(cComponent) = 1.4172D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.586988000D-10 ! Meter
          cWellDepth(cComponent) = 394.6076814D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.459306342D0
          cAspectRatio(cComponent) = 2.572671403D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.313467149D-10 ! Meter
          cWellDepth(cComponent) = 350.1835977D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.590582700D0
          cAspectRatio(cComponent) = 2.691729052D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.472582928D-10 ! Meter
            cWellDepth(cComponent) = 386.4513977D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.526709024D0
            cAspectRatio(cComponent) = 2.065532776D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 1.965610207D-10 ! Meter
                cWellDepth(cComponent) = 292.9580561D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.5349625157D0
                cAspectRatio(cComponent) = 3.884468640D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.719704124D-10 ! Meter
          cWellDepth(cComponent) = 479.0014247D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.448823527D0
          cAspectRatio(cComponent) = 2.853485380D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "CC4H6" .OR. cName == "CYCLOBUTENE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.748171234D-10 ! Meter
          cWellDepth(cComponent) = 483.9401652D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.433523426D0
          cAspectRatio(cComponent) = 2.103410466D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.897936424D-10 ! Meter
            cWellDepth(cComponent) = 497.6806784D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.420089948D0
            cAspectRatio(cComponent) = 1.648222929D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.192575733D-10 ! Meter
                cWellDepth(cComponent) = 532.1507846D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3164431381D0
                cAspectRatio(cComponent) = 1.193527136D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 54.0904D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.747578249D-10 ! Meter
          cWellDepth(cComponent) = 568.6772330D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.438552020D0
          cAspectRatio(cComponent) = 3.415421781D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "CC5H10" .OR. cName == "CYCLOPENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.788571988D-10 ! Meter
          cWellDepth(cComponent) = 579.9449111D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.437968397D0
          cAspectRatio(cComponent) = 2.716177713D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.982295524D-10 ! Meter
            cWellDepth(cComponent) = 604.8120409D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.414143081D0
            cAspectRatio(cComponent) = 2.044057794D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 3.517437373D-10 ! Meter
                cWellDepth(cComponent) = 673.6440308D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2685891167D0
                cAspectRatio(cComponent) = 1.301500826D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.451363260D-10 ! Meter
                cWellDepth(cComponent) = 671.0101972D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2814063385D0
                cAspectRatio(cComponent) = 1.295227762D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 70.1329D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.810093952D-10 ! Meter
          cWellDepth(cComponent) = 652.8130265D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.426632898D0
          cAspectRatio(cComponent) = 3.810777219D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "CC6H12" .OR. cName == "CYCLOHEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 2.846020235D-10 ! Meter
          cWellDepth(cComponent) = 611.6519091D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.469102602D0
          cAspectRatio(cComponent) = 2.989845561D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 3.069084160D-10 ! Meter
            cWellDepth(cComponent) = 657.0999700D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.425455924D0
            cAspectRatio(cComponent) = 2.217056061D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 3.700320900D-10 ! Meter
                cWellDepth(cComponent) = 752.8495621D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2592553010D0
                cAspectRatio(cComponent) = 1.342030199D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.610648054D-10 ! Meter
                cWellDepth(cComponent) = 758.2716704D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2687288599D0
                cAspectRatio(cComponent) = 1.360177656D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 84.1595D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.863690714D-10 ! Meter
          cWellDepth(cComponent) = 745.8373826D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.407294962D0
          cAspectRatio(cComponent) = 4.213862263D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.893606803D-10 ! Meter
          cWellDepth(cComponent) = 817.4183498D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.406295702D0
          cAspectRatio(cComponent) = 4.639115710D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.364724067D-10 ! Meter
          cWellDepth(cComponent) = 242.1438235D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456622914D0
          cAspectRatio(cComponent) = 2.406145387D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.917731224D-10 ! Meter
                cWellDepth(cComponent) = 284.5197866D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2841213690D0
                cAspectRatio(cComponent) = 1.220279190D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.489184869D-10 ! Meter
          cWellDepth(cComponent) = 339.9788780D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456218080D0
          cAspectRatio(cComponent) = 3.265836717D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.702983952D-10 ! Meter
          cWellDepth(cComponent) = 467.6032885D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.387548223D0
          cAspectRatio(cComponent) = 3.624922783D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN ! 15K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 1.765348503D-10 ! Meter
          cWellDepth(cComponent) = 313.3355532D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.543113051D0
          cAspectRatio(cComponent) = 3.604401315D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 1.898012536D-10 ! Meter
            cWellDepth(cComponent) = 355.1313076D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.468358858D0
            cAspectRatio(cComponent) = 2.722382686D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 1.674867829D-10 ! Meter
                cWellDepth(cComponent) = 334.1051578D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3118401261D0
                cAspectRatio(cComponent) = 4.333274466D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 1.587125140D-10 ! Meter
                cWellDepth(cComponent) = 339.7435404D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3439537109D0
                cAspectRatio(cComponent) = 4.339962871D0
              ELSE
                cDiameter(cComponent) = 1.714573890D-10 ! Meter
                cWellDepth(cComponent) = 351.0275683D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2836452641D0
                cAspectRatio(cComponent) = 4.127756559D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.645520833D-10 ! Meter
          cWellDepth(cComponent) = 617.0975616D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.478269449D0
          cAspectRatio(cComponent) = 3.072047352D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.857252284D-10 ! Meter
            cWellDepth(cComponent) = 663.5783367D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.433199057D0
            cAspectRatio(cComponent) = 2.267743227D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 3.443078169D-10 ! Meter
                cWellDepth(cComponent) = 763.1149114D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2616769254D0
                cAspectRatio(cComponent) = 1.383511002D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 3.350199214D-10 ! Meter
                cWellDepth(cComponent) = 769.2980334D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.2714405906D0
                cAspectRatio(cComponent) = 1.411989594D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.683264328D-10 ! Meter
          cWellDepth(cComponent) = 730.0764269D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.441918361D0
          cAspectRatio(cComponent) = 3.697187738D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.493580816D-10 ! Meter
          cWellDepth(cComponent) = 265.6127532D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.485946474D0
          cAspectRatio(cComponent) = 1.760108800D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.609234128D-10 ! Meter
            cWellDepth(cComponent) = 273.1643258D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.471847980D0
            cAspectRatio(cComponent) = 1.415449373D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.709954067D-10 ! Meter
                cWellDepth(cComponent) = 289.6830335D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3604302054D0
                cAspectRatio(cComponent) = 1.336613410D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.658982350D-10 ! Meter
                cWellDepth(cComponent) = 281.3543074D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3875803072D0
                cAspectRatio(cComponent) = 1.337674245D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.531420454D-10 ! Meter
          cWellDepth(cComponent) = 374.9806388D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.477593279D0
          cAspectRatio(cComponent) = 2.471139822D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.691934116D-10 ! Meter
            cWellDepth(cComponent) = 390.4138850D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.454300810D0
            cAspectRatio(cComponent) = 1.893898980D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          IF( PYHCBCorrectionLogical ) THEN ! Percus-Yevick correction
            IF( UseA1ForA2Logical ) THEN ! Use the same coefficient matrix of the effective packing fraction coefficients for both a1 and a2
              IF( ReferenceBoublikLogical ) THEN ! Replace the reference term with Boublik's expression
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                cDiameter(cComponent) = 2.970413648D-10 ! Meter
                cWellDepth(cComponent) = 436.9164870D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3030219613D0
                cAspectRatio(cComponent) = 1.501086938D0
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          ELSE
            IF( UseA1ForA2Logical ) THEN
              IF( ReferenceBoublikLogical ) THEN
                cDiameter(cComponent) = 2.897377260D-10 ! Meter
                cWellDepth(cComponent) = 430.3728045D0 * cBoltzmann ! Joule
                cPotentialRange(cComponent) = 0.3230885229D0
                cAspectRatio(cComponent) = 1.513026944D0
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            ELSE
              IF( ReferenceBoublikLogical ) THEN
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              ELSE
                WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
                CALL EXIT(  )
              END IF
            END IF
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN ! 0K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.069966865D-10 ! Meter
          cWellDepth(cComponent) = 18.14467890D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.689038273D0
          cAspectRatio(cComponent) = 2.860040066D-4
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 3.141542614D-10 ! Meter
          cWellDepth(cComponent) = 137.9121366D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.615774365D0
          cAspectRatio(cComponent) = 0.8861990317D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 3.091449175D-10 ! Meter
            cWellDepth(cComponent) = 137.7535712D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.636420795D0
            cAspectRatio(cComponent) = 0.8862225596D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN ! 10K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.751783756D-10 ! Meter
          cWellDepth(cComponent) = 288.4317154D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.495754517D0
          cAspectRatio(cComponent) = 1.992247275D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.755957822D-10 ! Meter
          cWellDepth(cComponent) = 394.6057319D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.459310920D0
          cAspectRatio(cComponent) = 2.679293508D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 5.296800976D-10 ! Meter
          cWellDepth(cComponent) = 350.1836236D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.590582710D0
          cAspectRatio(cComponent) = 0.2798213763D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 4.813117041D-10 ! Meter
            cWellDepth(cComponent) = 386.4529925D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.526706991D0
            cAspectRatio(cComponent) = 0.3704128340D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.875452121D-10 ! Meter
          cWellDepth(cComponent) = 479.0021228D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.448821698D0
          cAspectRatio(cComponent) = 2.978576269D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "CC4H6" .OR. cName == "CYCLOBUTENE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 5.410244305D-10 ! Meter
          cWellDepth(cComponent) = 483.9429531D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.433517032D0
          cAspectRatio(cComponent) = 0.3630580815D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 4.876906918D-10 ! Meter
            cWellDepth(cComponent) = 497.6806710D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.420090152D0
            cAspectRatio(cComponent) = 0.4856945301D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 54.0904D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.874028130D-10 ! Meter
          cWellDepth(cComponent) = 568.6770186D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.438552680D0
          cAspectRatio(cComponent) = 3.566639676D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "CC5H10" .OR. cName == "CYCLOPENTANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 6.419781578D-10 ! Meter
          cWellDepth(cComponent) = 579.9449089D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.437968227D0
          cAspectRatio(cComponent) = 0.2772470985D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 5.767653839D-10 ! Meter
            cWellDepth(cComponent) = 604.8121847D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.414142650D0
            cAspectRatio(cComponent) = 0.3747486996D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 70.1329D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.923854577D-10 ! Meter
          cWellDepth(cComponent) = 652.8141539D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.426628915D0
          cAspectRatio(cComponent) = 3.974913659D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "CC6H12" .OR. cName == "CYCLOHEXANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          cDiameter(cComponent) = 6.945876171D-10 ! Meter
          cWellDepth(cComponent) = 611.6508954D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.469104640D0
          cAspectRatio(cComponent) = 0.2515352663D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 6.241624186D-10 ! Meter
            cWellDepth(cComponent) = 657.0995811D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.425456711D0
            cAspectRatio(cComponent) = 0.3428360799D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 84.1595D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.966951621D-10 ! Meter
          cWellDepth(cComponent) = 745.8377800D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.407294172D0
          cAspectRatio(cComponent) = 4.388483057D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.987201893D-10 ! Meter
          cWellDepth(cComponent) = 817.4231533D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.406291020D0
          cAspectRatio(cComponent) = 4.822539398D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.533396220D-10 ! Meter
          cWellDepth(cComponent) = 242.1442215D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456620830D0
          cAspectRatio(cComponent) = 2.499017073D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN ! 20K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.610036528D-10 ! Meter
          cWellDepth(cComponent) = 339.9787447D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.456218378D0
          cAspectRatio(cComponent) = 3.411148957D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.818962471D-10 ! Meter
          cWellDepth(cComponent) = 467.6031656D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.387547234D0
          cAspectRatio(cComponent) = 3.783397022D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN ! 15K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 1.841593007D-10 ! Meter
          cWellDepth(cComponent) = 313.3355137D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.543113252D0
          cAspectRatio(cComponent) = 3.762242139D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.013320685D-10 ! Meter
            cWellDepth(cComponent) = 355.1313175D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.468358992D0
            cAspectRatio(cComponent) = 2.839462306D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(2G0)" ) "There are two set of parameters for BENZENE. Set [1] has an aspect ratio of a (DISK) and ", &
          &                   "set [2] has an aspect ratio of a (ROD). Select the dataset for BENZENE..."
          READ( *, * ) DataSet
          IF( DataSet == 1 ) THEN
            cDiameter(cComponent) = 6.563947943D-10 ! Meter
            cWellDepth(cComponent) = 617.0966870D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.478271203D0
            cAspectRatio(cComponent) = 0.2447707866D0
            IF( ReferenceBoublikLogical ) THEN
              cDiameter(cComponent) = 5.892106083D-10 ! Meter
              cWellDepth(cComponent) = 663.5781624D0 * cBoltzmann ! Joule
              cPotentialRange(cComponent) = 1.433199333D0
              cAspectRatio(cComponent) = 0.3346221628D0
            END IF
          ELSE IF( DataSet == 2 ) THEN
            cDiameter(cComponent) = 2.783820373D-10 ! Meter
            cWellDepth(cComponent) = 617.0964658D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.478271676D0
            cAspectRatio(cComponent) = 3.208714557D0
          ELSE
            WRITE( *, "(G0)" ) "Invalid dataset. Exiting..."
            CALL EXIT(  )
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN ! 40K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(2G0)" ) "There are two set of parameters for TOLUENE. Set [1] has an aspect ratio of 0.245 (DISK) and ", &
          &                   "set [2] has an aspect ratio of 3.209 (ROD). Select the dataset for TOLUENE..."
          READ( *, * ) DataSet
          IF( DataSet == 1 ) THEN
            cDiameter(cComponent) = 7.453797472D-10 ! Meter
            cWellDepth(cComponent) = 730.0760905D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441919089D0
            cAspectRatio(cComponent) = 0.2035762393D0
          ELSE IF( DataSet == 2 ) THEN
            cDiameter(cComponent) = 2.795754772D-10 ! Meter
            cWellDepth(cComponent) = 730.0761976D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441918977D0
            cAspectRatio(cComponent) = 3.858004931D0
          ELSE
            WRITE( *, "(G0)" ) "Invalid dataset. Exiting..."
            CALL EXIT(  )
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(2G0)" ) "There are two set of parameters for TOLUENE. Set [1] has an aspect ratio of 0.245 (DISK) and ", &
          &                   "set [2] has an aspect ratio of 3.209 (ROD). Select the dataset for TOLUENE..."
          READ( *, * ) DataSet
          IF( DataSet == 1 ) THEN
            cDiameter(cComponent) = 7.453797472D-10 ! Meter
            cWellDepth(cComponent) = 730.0760905D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441919089D0
            cAspectRatio(cComponent) = 0.2035762393D0
          ELSE IF( DataSet == 2 ) THEN
            cDiameter(cComponent) = 2.795754772D-10 ! Meter
            cWellDepth(cComponent) = 730.0761976D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.441918977D0
            cAspectRatio(cComponent) = 3.858004931D0
          ELSE
            WRITE( *, "(G0)" ) "Invalid dataset. Exiting..."
            CALL EXIT(  )
          END IF
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.775749233D-10 ! Meter
          cWellDepth(cComponent) = 265.6128236D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.485946170D0
          cAspectRatio(cComponent) = 1.759379077D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 3.072062853D-10 ! Meter
            cWellDepth(cComponent) = 273.1664620D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.471837703D0
            cAspectRatio(cComponent) = 1.275734754D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.705693739D-10 ! Meter
          cWellDepth(cComponent) = 374.9800261D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.477595204D0
          cAspectRatio(cComponent) = 2.569697576D0
          IF( ReferenceBoublikLogical ) THEN
            cDiameter(cComponent) = 2.962458126D-10 ! Meter
            cWellDepth(cComponent) = 390.4139753D0 * cBoltzmann ! Joule
            cPotentialRange(cComponent) = 1.454300504D0
            cAspectRatio(cComponent) = 1.921200537D0
          END IF
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN ! 0K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          cDiameter(cComponent) = 2.729416744D-10 ! Meter
          cWellDepth(cComponent) = 17.33438150D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.773550224D0
          cAspectRatio(cComponent) = 0.9053531042D0
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          cDiameter(cComponent) = 2.729416744D-10 ! Meter
          cWellDepth(cComponent) = 17.33438150D0 * cBoltzmann ! Joule
          cPotentialRange(cComponent) = 1.773550224D0
          cAspectRatio(cComponent) = 0.9053531042D0
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    END IF
  ELSE ! No Zhang correction
    IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-Revolution
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
      IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 16.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 30.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 44.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN ! 30K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 58.0782D-3 ! Kilograms / mol
      ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 72.0940D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 86.1096D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 100.1252D-3 ! Kilograms / mol
      ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 114.1408D-3 ! Kilograms / mol
      ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 87.9936D-3 ! Kilograms / mol
      ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 137.9904D-3 ! Kilograms / mol
      ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 187.9872D-3 ! Kilograms / mol
      ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 43.98983D-3 ! Kilograms / mol
      ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 78.0470D-3 ! Kilograms / mol
      ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 92.0626D-3 ! Kilograms / mol
      ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 28.0313D-3 ! Kilograms / mol
      ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN ! 25K below the critical temperature
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 42.0797D-3 ! Kilograms / mol
      ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 2.01568D-3 ! Kilograms / mol
      ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 52.0238D-3 ! Kilograms / mol
      ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
        IF( PotentialTypeLogical(1) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(2) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(3) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        ELSE IF( PotentialTypeLogical(4) ) THEN
          WRITE( *, "(G0)" ) "Properties not available for the selected potential. Exiting..."
          CALL EXIT(  )
        END IF
        cLength(cComponent) = cAspectRatio(cComponent) * cDiameter(cComponent)
        cMolarMass(cComponent) = 120.0220D-3 ! Kilograms / mol
      ELSE
        WRITE( *, "(G0)" ) "Selected molecule does not exist in the database. Exiting..."
        CALL EXIT(  )
      END IF
    END IF
  END IF
ELSE ! Including higher-order TPT terms
  WRITE( *, "(G0)" ) "Properties not available for higher order terms. Exiting..."
  CALL EXIT(  )
END IF

! Diameter of the equivalent sphere
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( 1.D0 + 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
END IF

! Field properties (CSW)
fDiameter(cComponent) = cDiameter(cComponent) + cPotentialRange(cComponent) * cDiameterSphere(cComponent)
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  fLength(cComponent) = cLength(cComponent) + cPotentialRange(cComponent) * cDiameterSphere(cComponent)
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  fLength(cComponent) = cLength(cComponent)
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  fLength(cComponent) = cLength(cComponent)
END IF
fAspectRatio(cComponent) = fLength(cComponent) / fDiameter(cComponent)

RETURN

END SUBROUTINE Properties

! ************************************************************************************************ !
!                Defines all specific heat parameters of molecules within database                 !
! ************************************************************************************************ !
SUBROUTINE Specific_Heat_Parameters( cComponent, cName )

! Use some variables of the global variables module
USE GlobalVariables, ONLY: Int64, Real64, cBoltzmann, SpecificHeatReference, TemperatureParameter

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( LEN= * ) :: cName ! Component name

! Thermodynamics Research Center Polynomial Equation: Temperature range (50-1000 K)
IF( SpecificHeatReference == "TRC" ) THEN
  IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
    TemperatureParameter(cComponent,1) = 4.568D0
    TemperatureParameter(cComponent,2) = -8.975D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 3.631D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -3.407D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 1.091D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 4.178D0
    TemperatureParameter(cComponent,2) = -4.427D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 5.660D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -6.651D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 2.487D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN
    TemperatureParameter(cComponent,1) = 3.847D0
    TemperatureParameter(cComponent,2) = 5.131D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 6.011D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -7.893D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 3.079D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN
    TemperatureParameter(cComponent,1) = 4.493D0
    TemperatureParameter(cComponent,2) = 18.097D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 12.744D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = 16.049D-8  ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 6.426D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN
    TemperatureParameter(cComponent,1) = 5.547D0
    TemperatureParameter(cComponent,2) = 5.536D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 8.057D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -10.571D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 4.134D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CC4H6" .OR. cName == "CYCLOBUTENE" ) THEN
    WRITE( *, "(2G0)" ) "Parameters of the Thermodynamics Research Center (TRC) Polynomial Equation for selected molecule are ", &
    &                   "not available. Using parameters for BUTANE..."
    TemperatureParameter(cComponent,1) = 5.547D0
    TemperatureParameter(cComponent,2) = 5.536D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 8.057D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -10.571D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 4.134D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN
    TemperatureParameter(cComponent,1) = 7.554D0
    TemperatureParameter(cComponent,2) = -0.368D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 11.846D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -14.939D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 5.753D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CC5H10" .OR. cName == "CYCLOPENTANE" ) THEN
    TemperatureParameter(cComponent,1) = 5.019D0
    TemperatureParameter(cComponent,2) = -19.734D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 17.917D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -21.696D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 8.215D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN
    TemperatureParameter(cComponent,1) = 8.831D0
    TemperatureParameter(cComponent,2) = -0.166D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 14.302D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -18.3134D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 7.124D-11   ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CC6H12" .OR. cName == "CYCLOHEXANE" ) THEN
    TemperatureParameter(cComponent,1) = 4.035D0
    TemperatureParameter(cComponent,2) = -4.433D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 16.834D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -20.775D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 7.746D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN
    TemperatureParameter(cComponent,1) = 9.634D0
    TemperatureParameter(cComponent,2) = 4.156D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 15.494D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -20.066D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 7.770D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN
    TemperatureParameter(cComponent,1) = 10.824D0
    TemperatureParameter(cComponent,2) = 4.883D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 17.751D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -23.137D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 8.980D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 2.643D0
    TemperatureParameter(cComponent,2) = 15.383D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 0.850D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -2.940D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 1.469D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 2.525D0
    TemperatureParameter(cComponent,2) = 43.543D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = -2.948D-5 ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -0.630D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 0.967D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN
    TemperatureParameter(cComponent,1) = 1.605D0
    TemperatureParameter(cComponent,2) = 76.488D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = -8.707D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = 4.540D-8   ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = -0.856D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
    TemperatureParameter(cComponent,1) = 3.259D0
    TemperatureParameter(cComponent,2) = 1.356D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 1.502D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -2.374D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 1.056D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN
    TemperatureParameter(cComponent,1) = 3.551D0
    TemperatureParameter(cComponent,2) = -6.184D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 14.365D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -19.807D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 8.234D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN
    TemperatureParameter(cComponent,1) = 3.866D0
    TemperatureParameter(cComponent,2) = 3.558D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 13.356D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -18.659D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 7.690D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
    TemperatureParameter(cComponent,1) = 4.221D0
    TemperatureParameter(cComponent,2) = -8.782D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 5.795D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -6.729D-8  ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 2.511D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN
    TemperatureParameter(cComponent,1) = 3.834D0
    TemperatureParameter(cComponent,2) = 3.893D-3  ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 4.688D-5   ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = 6.013D-8  ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 2.283D-11  ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
    TemperatureParameter(cComponent,1) = 2.883D0
    TemperatureParameter(cComponent,2) = 3.681D-3   ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = -0.772D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = 0.692D-8   ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = -0.213D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 4.150D0
    TemperatureParameter(cComponent,2) = -5.584D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = 4.384D-5  ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -5.160D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 1.920D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 3.146D0
    TemperatureParameter(cComponent,2) = 29.937D-3 ! 1 / Kelvin
    TemperatureParameter(cComponent,3) = -0.056D-5 ! 1 / Kelvin²
    TemperatureParameter(cComponent,4) = -3.019D-8 ! 1 / Kelvin³
    TemperatureParameter(cComponent,5) = 1.669D-11 ! 1 / (Kelvin² . Kelvin²)
  ELSE
    WRITE( *, "(2G0)" ) "Parameters of the Thermodynamics Research Center (TRC) Polynomial Equation for selected molecule are ", &
    &                   "not available. Exiting..."
    CALL EXIT(  )
  END IF
! Shomate equation: Temperature range (298-1300 K)
ELSE IF( SpecificHeatReference == "NIST" ) THEN
  IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
    TemperatureParameter(cComponent,1) = -0.703029D0 ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 108.4773D0  ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -42.52157D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 5.862788D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = 0.678565D0  ! Joule . Kelvin / mol
  ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 15.96778D0  ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 210.3318D0  ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -189.4657D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 62.20227D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = -0.217317D0 ! Joule . Kelvin / mol
  ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
    TemperatureParameter(cComponent,1) = 69.74728D0  ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 225.3599D0  ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -174.1628D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 47.14657D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = -1.455677D0 ! Joule . Kelvin / mol
  ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
    TemperatureParameter(cComponent,1) = 24.99735D0  ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 55.18696D0  ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -33.69137D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 7.948387D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = -0.136638D0 ! Joule . Kelvin / mol
  ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
    TemperatureParameter(cComponent,1) = -6.387880D0 ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 184.4019D0  ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -112.9718D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 28.49593D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = 0.315540D0  ! Joule . Kelvin / mol
  ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
    TemperatureParameter(cComponent,1) = 33.066178D0  ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = -11.363417D0 ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = 11.432816D0  ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = -2.772874D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = -0.158558D0  ! Joule . Kelvin / mol
  ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
    TemperatureParameter(cComponent,1) = -6.098682D0 ! Joule / (mol . Kelvin)
    TemperatureParameter(cComponent,2) = 179.22D0    ! Joule / (mol . Kelvin²)
    TemperatureParameter(cComponent,3) = -122.3682D0 ! Joule / (mol . Kelvin³)
    TemperatureParameter(cComponent,4) = 32.30207D0  ! Joule / (mol . Kelvin⁴)
    TemperatureParameter(cComponent,5) = 0.491361D0  ! Joule . Kelvin / mol
  ELSE
    WRITE( *, "(G0)" ) "Shomate equation parameters for selected molecule are not available. Exiting..."
    CALL EXIT(  )
  END IF
END IF

RETURN

END SUBROUTINE Specific_Heat_Parameters

! ************************************************************************************************ !
!                Defines all specific heat parameters of molecules within database                 !
! ************************************************************************************************ !
SUBROUTINE Structural_Formula( cComponent, cName )

! Use some variables of the global variables module
USE GlobalVariables, ONLY: Int64, Real64, cFormulaName

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( LEN= * ) :: cName ! Component name

! Formulas
IF( cName == "CH4" .OR. cName == "METHANE" ) THEN
  cFormulaName(cComponent) = "CH4"
ELSE IF( cName == "C2H6" .OR. cName == "ETHANE" ) THEN
  cFormulaName(cComponent) = "C2H6"
ELSE IF( cName == "C3H8" .OR. cName == "PROPANE" ) THEN
  cFormulaName(cComponent) = "C3H8"
ELSE IF( cName == "CC3H6" .OR. cName == "CYCLOPROPANE" ) THEN
  cFormulaName(cComponent) = "CC3H6"
ELSE IF( cName == "C4H10" .OR. cName == "BUTANE" ) THEN
  cFormulaName(cComponent) = "C4H10"
ELSE IF( cName == "CC4H6" .OR. cName == "CYCLOBUTENE" ) THEN
  cFormulaName(cComponent) = "CC4H6"
ELSE IF( cName == "C5H12" .OR. cName == "PENTANE" ) THEN
  cFormulaName(cComponent) = "C5H12"
ELSE IF( cName == "CC5H10" .OR. cName == "CYCLOPENTANE" ) THEN
  cFormulaName(cComponent) = "CC5H10"
ELSE IF( cName == "C6H14" .OR. cName == "HEXANE" ) THEN
  cFormulaName(cComponent) = "C6H14"
ELSE IF( cName == "CC6H12" .OR. cName == "CYCLOHEXANE" ) THEN
  cFormulaName(cComponent) = "CC6H12"
ELSE IF( cName == "C7H16" .OR. cName == "HEPTANE" ) THEN
  cFormulaName(cComponent) = "C7H16"
ELSE IF( cName == "C8H18" .OR. cName == "OCTANE" ) THEN
  cFormulaName(cComponent) = "C8H18"
ELSE IF( cName == "CF4" .OR. cName == "TETRAFLUOROMETHANE" ) THEN
  cFormulaName(cComponent) = "CF4"
ELSE IF( cName == "C2F6" .OR. cName == "HEXAFLUOROETHANE" ) THEN
  cFormulaName(cComponent) = "C2F6"
ELSE IF( cName == "C3F8" .OR. cName == "OCTAFLUOROPROPANE" ) THEN
  cFormulaName(cComponent) = "C3F8"
ELSE IF( cName == "CO2" .OR. cName == "CARBONDIOXIDE" ) THEN
  cFormulaName(cComponent) = "CO2"
ELSE IF( cName == "C6H6" .OR. cName == "BENZENE" ) THEN
  cFormulaName(cComponent) = "C6H6"
ELSE IF( cName == "C7H8" .OR. cName == "TOLUENE" ) THEN
  cFormulaName(cComponent) = "C7H8"
ELSE IF( cName == "C2H4" .OR. cName == "ETHYLENE" ) THEN
  cFormulaName(cComponent) = "C2H4"
ELSE IF( cName == "C3H6" .OR. cName == "PROPENE" ) THEN
  cFormulaName(cComponent) = "C3H6"
ELSE IF( cName == "H2" .OR. cName == "HYDROGEN" ) THEN
  cFormulaName(cComponent) = "H2"
ELSE IF( cName == "R32" .OR. cName == "DIFLUOROMETHANE" ) THEN
  cFormulaName(cComponent) = "R32"
ELSE IF( cName == "R125" .OR. cName == "PENTAFLUOROETHANE" ) THEN
  cFormulaName(cComponent) = "R125"
END IF


RETURN

END SUBROUTINE Structural_Formula

! ************************************************************************************************ !
!                Defines all reduced/multiple molecular properties within database                 !
! ************************************************************************************************ !
SUBROUTINE Combined_Properties(  )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! ************************************************************************************************ !
! REDUCED PROPERTIES                                                                               !
! ************************************************************************************************ !

! Diameter
aDiameter = cDiameter / 1.D-10 ! Å
aDiameterField = fDiameter / 1.D-10 ! Å
ijaDiameter = ijDiameter / 1.D-10 ! Å

! Spherical diameter
aDiameterSphere = cDiameterSphere / 1.D-10 ! Å
ijaDiameterSphere = ijDiameterSphere / 1.D-10 ! Å
ijDiameterSphereCubic = ijDiameterSphere * ijDiameterSphere * ijDiameterSphere ! m³
ijaDiameterSphereCubic = ijaDiameterSphere * ijaDiameterSphere * ijaDiameterSphere ! Å³

! Length
aLength = cLength / 1.D-10 ! Å
aLengthField = fLength / 1.D-10 ! Å
ijaLength = ijLength / 1.D-10 ! Å

! Potential depth
aWellDepth = cWellDepth / cBoltzmann ! K
ijaWellDepth = ijWellDepth / cBoltzmann ! K

! Potential range
ijPotentialRangeSquared = ijPotentialRange * ijPotentialRange
ijPotentialRangeCubic = ijPotentialRange * ijPotentialRange * ijPotentialRange

RETURN

END SUBROUTINE Combined_Properties

! ************************************************************************************************ !
! Isihara-Hadwiger Theorem (Second Virial Coefficient)                                             !
! ************************************************************************************************ !
SUBROUTINE Isihara_Hadwiger_Theorem(  )

! Uses two modules: global variables
USE GlobalVariables

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Component indexes

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: MolecularVolume  ! Å³
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: fMolecularVolume ! Å³
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: SurfaceArea      ! Å²
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: fSurfaceArea     ! Å²
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: CurvatureRadius  ! Å
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: fCurvatureRadius ! Å
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: Radius           ! Å
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: fRadius          ! Å

! Radius of hard convex body [Å]
Radius = 0.5D0 * aDiameter

! Isihara-Hadwiger theorem
DO iComponent = 1, nComponents
  IF( GeometrySpecification(iComponent,1) ) THEN ! Ellipsoids-of-revolution
    ! Molecular volume [Å³]
    MolecularVolume(iComponent) = ( cPi / 6.D0 ) * aDiameter(iComponent) * aDiameter(iComponent) * aDiameter(iComponent) * &
    &     cAspectRatio(iComponent)
    ! Prolate ellipsoids-of-revolution
    IF( cAspectRatio(iComponent) > 1.D0 ) THEN
      ! Radius of curvature [Å]
      CurvatureRadius(iComponent) = 0.5D0 * cAspectRatio(iComponent) * Radius(iComponent) + ( 0.5D0 * Radius(iComponent) / &
      &     DSQRT( cAspectRatio(iComponent) * cAspectRatio(iComponent) - 1.D0 ) ) * DATANH( DSQRT( cAspectRatio(iComponent) * &
      &     cAspectRatio(iComponent) - 1.D0 ) / cAspectRatio(iComponent) )
      ! Surface area [Å²]
      SurfaceArea(iComponent) = ( 2.D0 * cPi * Radius(iComponent) * Radius(iComponent) / DSQRT( cAspectRatio(iComponent) * &
      &     cAspectRatio(iComponent) - 1.D0 ) ) * ( ( cAspectRatio(iComponent) * cAspectRatio(iComponent) * DASIN( DSQRT( &
      &     cAspectRatio(iComponent) * cAspectRatio(iComponent) - 1.D0 ) / cAspectRatio(iComponent) ) ) + ( DSQRT( &
      &     cAspectRatio(iComponent) * cAspectRatio(iComponent) - 1.D0 ) ) )
    ! Oblate ellipsoids-of-revolution
    ELSE IF( cAspectRatio(iComponent) < 1.D0 ) THEN
      ! Radius of curvature [Å]
      CurvatureRadius(iComponent) = 0.5D0 * cAspectRatio(iComponent) * Radius(iComponent) + ( 0.5D0 * Radius(iComponent) / &
      &     DSQRT( 1.D0 - cAspectRatio(iComponent) * cAspectRatio(iComponent) ) ) * DATAN( DSQRT( 1.D0 - cAspectRatio(iComponent) &
      &     * cAspectRatio(iComponent) ) / cAspectRatio(iComponent) )
      ! Surface area [Å²]
      SurfaceArea(iComponent) = ( 2.D0 * cPi * Radius(iComponent) * Radius(iComponent) / DSQRT( 1.D0 - cAspectRatio(iComponent) * &
      &     cAspectRatio(iComponent) ) ) * ( ( cAspectRatio(iComponent) * cAspectRatio(iComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 - &
      &     cAspectRatio(iComponent) * cAspectRatio(iComponent) ) ) / cAspectRatio(iComponent) ) ) + ( DSQRT( 1.D0 - &
      &     cAspectRatio(iComponent) * cAspectRatio(iComponent) ) ) )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      CurvatureRadius(iComponent) = Radius(iComponent)
      ! Surface area [Å²]
      SurfaceArea(iComponent) = 4.D0 * cPi * Radius(iComponent) * Radius(iComponent)
    END IF
  ELSE IF( GeometrySpecification(iComponent,2) ) THEN ! Spherocylinders
    ! Molecular volume [Å³]
    MolecularVolume(iComponent) = ( cPi / 4.D0 ) * aDiameter(iComponent) * aDiameter(iComponent) * aDiameter(iComponent) * &
    &     ( cAspectRatio(iComponent) + (2.D0 / 3.D0) )
    ! Prolate spherocylinders
    IF( cAspectRatio(iComponent) > 0.D0 ) THEN
      ! Radius of curvature [Å]
      CurvatureRadius(iComponent) = 0.5D0 * ( cAspectRatio(iComponent) + 2.D0 ) * Radius(iComponent)
      ! Surface area [Å²]
      SurfaceArea(iComponent) = 4.D0 * cPi * Radius(iComponent) * Radius(iComponent) * ( cAspectRatio(iComponent) + 1.D0 )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      CurvatureRadius(iComponent) = Radius(iComponent)
      ! Surface area [Å²]
      SurfaceArea(iComponent) = 4.D0 * cPi * Radius(iComponent) * Radius(iComponent)
    END IF
  ELSE IF( GeometrySpecification(iComponent,3) ) THEN ! Cylinders
    ! Molecular volume [Å³]
    MolecularVolume(iComponent) = ( cPi / 4.D0 ) * aDiameter(iComponent) * aDiameter(iComponent) * aDiameter(iComponent) * &
    &     cAspectRatio(iComponent)
    ! Radius of curvature [Å]
    CurvatureRadius(iComponent) = 0.5D0 * ( cAspectRatio(iComponent) + (0.5D0 * cPi) ) * Radius(iComponent)
    ! Surface area [Å²]
    SurfaceArea(iComponent) = 2.D0 * cPi * Radius(iComponent) * Radius(iComponent) * ( 1.D0 + 2.D0 * cAspectRatio(iComponent) )
  END IF
END DO

! Radius
fRadius = 0.5D0 * aDiameterField

! Morphological descriptors (field)
DO iComponent = 1, nComponents
  IF( GeometrySpecification(iComponent,1) ) THEN ! Ellipsoids-of-revolution
    ! Molecular volume [Å³]
    fMolecularVolume(iComponent) = ( cPi / 6.D0 ) * aDiameterField(iComponent) * aDiameterField(iComponent) * &
    &     aDiameterField(iComponent) * fAspectRatio(iComponent)
    ! Prolate ellipsoids-of-revolution
    IF( fAspectRatio(iComponent) > 1.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(iComponent) = 0.5D0 * fAspectRatio(iComponent) * fRadius(iComponent) + ( 0.5D0 * fRadius(iComponent) / &
      &     DSQRT( fAspectRatio(iComponent) * fAspectRatio(iComponent) - 1.D0 ) ) * DATANH( DSQRT( fAspectRatio(iComponent) * &
      &     fAspectRatio(iComponent) - 1.D0 ) / fAspectRatio(iComponent) )
      ! Surface area [Å²]
      fSurfaceArea(iComponent) = ( 2.D0 * cPi * fRadius(iComponent) * fRadius(iComponent) / DSQRT( fAspectRatio(iComponent) * &
      &     fAspectRatio(iComponent) - 1.D0 ) ) * ( ( fAspectRatio(iComponent) * fAspectRatio(iComponent) * DASIN( DSQRT( &
      &     fAspectRatio(iComponent) * fAspectRatio(iComponent) - 1.D0 ) / fAspectRatio(iComponent) ) ) + ( DSQRT( &
      &     fAspectRatio(iComponent) * fAspectRatio(iComponent) - 1.D0 ) ) )
    ! Oblate ellipsoids-of-revolution
    ELSE IF( fAspectRatio(iComponent) < 1.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(iComponent) = 0.5D0 * fAspectRatio(iComponent) * fRadius(iComponent) + ( 0.5D0 * fRadius(iComponent) / &
      &     DSQRT( 1.D0 - fAspectRatio(iComponent) * fAspectRatio(iComponent) ) ) * DATAN( DSQRT( 1.D0 - fAspectRatio(iComponent) &
      &     * fAspectRatio(iComponent) ) / fAspectRatio(iComponent) )
      ! Surface area [Å²]
      fSurfaceArea(iComponent) = ( 2.D0 * cPi * fRadius(iComponent) * fRadius(iComponent) / DSQRT( 1.D0 - fAspectRatio(iComponent) &
      &     * fAspectRatio(iComponent) ) ) * ( ( fAspectRatio(iComponent) * fAspectRatio(iComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 &
      &     - fAspectRatio(iComponent) * fAspectRatio(iComponent) ) ) / fAspectRatio(iComponent) ) ) + ( DSQRT( 1.D0 - &
      &     fAspectRatio(iComponent) * fAspectRatio(iComponent) ) ) )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      fCurvatureRadius(iComponent) = fRadius(iComponent)
      ! Surface area [Å²]
      fSurfaceArea(iComponent) = 4.D0 * cPi * fRadius(iComponent) * fRadius(iComponent)
    END IF
  ELSE IF( GeometrySpecification(iComponent,2) .OR. GeometrySpecification(iComponent,3) ) THEN ! Spherocylinders and cylinders
    ! Molecular volume [Å³]
    fMolecularVolume(iComponent) = ( cPi / 4.D0 ) * aDiameterField(iComponent) * aDiameterField(iComponent) * &
    &     aDiameterField(iComponent) * ( fAspectRatio(iComponent) + (2.D0 / 3.D0) )
    ! Prolate spherocylinders
    IF( fAspectRatio(iComponent) > 0.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(iComponent) = 0.5D0 * ( fAspectRatio(iComponent) + 2.D0 ) * fRadius(iComponent)
      ! Surface area [Å²]
      fSurfaceArea(iComponent) = 4.D0 * cPi * fRadius(iComponent) * fRadius(iComponent) * ( fAspectRatio(iComponent) + 1.D0 )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      fCurvatureRadius(iComponent) = fRadius(iComponent)
      ! Surface area [Å²]
      fSurfaceArea(iComponent) = 4.D0 * cPi * fRadius(iComponent) * fRadius(iComponent)
    END IF
  END IF
END DO

! Isihara-Hadwiger second virial coefficients
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijSecondVirialCoefficient(iComponent,jComponent) = 0.5D0 * ( MolecularVolume(iComponent) + MolecularVolume(jComponent) ) + &
    &     0.5D0 * ( CurvatureRadius(iComponent) * SurfaceArea(jComponent) ) + 0.5D0 * ( CurvatureRadius(jComponent) * &
    &     SurfaceArea(iComponent) )
    ijSecondVirialCoefficientField(iComponent,jComponent) = 0.5D0 * ( fMolecularVolume(iComponent) + &
    &     fMolecularVolume(jComponent) ) + 0.5D0 * ( fCurvatureRadius(iComponent) * fSurfaceArea(jComponent) ) + 0.5D0 * &
    &     ( fCurvatureRadius(jComponent) * fSurfaceArea(iComponent) )
    ijNonSphericity(iComponent,jComponent) = ( 6.D0 * ijSecondVirialCoefficient(iComponent,jComponent) / cPi / &
    &     ijaDiameterSphereCubic(iComponent,jComponent) - 1.D0 ) / 3.D0
  END DO
END DO

! Isihara-Hadwiger theorem (hard spheres)
DO iComponent = 1, nComponents
  ! Molecular volume [Å³]
  MolecularVolume(iComponent) = ( cPi / 6.D0 ) * aDiameterSphere(iComponent) * aDiameterSphere(iComponent) * &
  &     aDiameterSphere(iComponent)
  ! Radius of curvature [Å]
  CurvatureRadius(iComponent) = 0.5D0 * aDiameterSphere(iComponent)
  ! Surface area [Å²]
  SurfaceArea(iComponent) = cPi * aDiameterSphere(iComponent) * aDiameterSphere(iComponent)
END DO

! Isihara-Hadwiger second virial coefficients (hard spheres)
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijHSSecondVirialCoefficient(iComponent,jComponent) = 0.5D0 * ( MolecularVolume(iComponent) + MolecularVolume(jComponent) ) + &
    &     0.5D0 * ( CurvatureRadius(iComponent) * SurfaceArea(jComponent) ) + 0.5D0 * ( CurvatureRadius(jComponent) * &
    &     SurfaceArea(iComponent) )
  END DO
END DO

! Ratio of the second virial coefficients (nonspherical and spherical)
ijRatioSecondVirialCoefficient = ijSecondVirialCoefficient / ijHSSecondVirialCoefficient

RETURN

END SUBROUTINE Isihara_Hadwiger_Theorem

END MODULE Substances
