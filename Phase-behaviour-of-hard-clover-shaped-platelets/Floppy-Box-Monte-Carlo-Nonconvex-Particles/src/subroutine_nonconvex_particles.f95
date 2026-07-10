! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!                   This code contains all subroutines used in the main program                   !
!                               to design the nonconvex particles.                                !
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

! *********************************************************************************************** !
!         This subroutine allows the user to choose the surface geometry of the molecules         !
! *********************************************************************************************** !
SUBROUTINE GeometrySelection(  )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
!  CLOVER-LIKE GEOMETRIES                                                                         !
!                                                                                                 !
!  The geometry of the cross sections consists of four circles and are constructed as follows:    !
!                                                                                                 !
!  (1) = RULES: A circle intersects the centers of the two nearest neighboring circles            !
!               A circle is equidistant from the two nearest neighboring circles                  !
!                                                                                                 !
!  (2) = RULES: All four circles intersect the center of mass of the composed geometry            !
!               A circle is equidistant from the two nearest neighboring circles                  !
! *********************************************************************************************** !
PetalArrangementLogical = .FALSE.

! *********************************************************************************************** !
! MOLECULAR PROPERTIES                                                                            !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )

! Skip
READ( 100, * ) Dummy, Dummy
! Cross-section geometry
READ( 100, * ) Dummy, Arrangement
! Rule (1)
IF( Arrangement == 1 ) THEN
  PetalArrangementLogical(1) = .TRUE.
! Rule (2)
ELSE IF( Arrangement == 2 ) THEN
  PetalArrangementLogical(2) = .TRUE.
END IF

CLOSE( 100 )

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
SUBROUTINE CrossSectionArea( SurfaceArea )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: AreaCircle                                         ! Area of the circles
REAL( Kind= Real64 ) :: AreaIntersection2Circles, AreaIntersection4Circles ! Area of the circular intersections
REAL( Kind= Real64 ) :: SurfaceArea                                        ! Surface area

! Arrangement[1]
IF( PetalArrangementLogical(1) ) THEN
  ! Area of circles
  AreaCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  AreaIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( 2.D0 * cPi / 3.D0 ) - ( DSQRT( 3.D0 ) / 2.D0 ) )
  ! Area of the intersection between four circles
  AreaIntersection4Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( cPi / 3.D0 ) - ( DSQRT( 3.D0 ) ) + 1.D0 )
  ! Surface area
  SurfaceArea = 4.D0 * AreaCircle - 4.D0 * AreaIntersection2Circles + AreaIntersection4Circles
! Arrangement[2]
ELSE IF( PetalArrangementLogical(2) ) THEN
  ! Area of circles
  AreaCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  AreaIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( cPi / 2.D0 - 1.D0 )
  ! Surface area
  SurfaceArea = 4.D0 * AreaCircle - 4.D0 * AreaIntersection2Circles
END IF

RETURN

END SUBROUTINE CrossSectionArea