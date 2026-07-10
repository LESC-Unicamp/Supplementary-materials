! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!    This code contains the subroutine used in the main program to create the list of images.     !
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
!                    This subroutine computes the images of a central particle                    !
! *********************************************************************************************** !
SUBROUTINE ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: xLayer, yLayer, zLayer ! X-, Y-, and Z- coordinates
INTEGER( Kind= Int64 ) :: cLayer                 ! Counter (layers)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 )          :: ScalingDistanceUnitBox      ! Scaled position, particles (unit box)
REAL( Kind= Real64 ), DIMENSION( 3, nImages ) :: ScalingDistanceImageUnitBox ! Scaled position, images (unit box)

! Position of images
cLayer = 1
DO zLayer = - nLayers, nLayers
  DO yLayer = - nLayers, nLayers
    DO xLayer = - nLayers, nLayers
      IF( zLayer == 0 .AND. yLayer == 0 .AND. xLayer == 0 ) CYCLE ! This is not an image
      ScalingDistanceImageUnitBox(1,cLayer) = ScalingDistanceUnitBox(1) + DBLE( xLayer )
      ScalingDistanceImageUnitBox(2,cLayer) = ScalingDistanceUnitBox(2) + DBLE( yLayer )
      ScalingDistanceImageUnitBox(3,cLayer) = ScalingDistanceUnitBox(3) + DBLE( zLayer )
      cLayer = cLayer + 1
    END DO
  END DO
END DO

RETURN

END SUBROUTINE ImageConstruction