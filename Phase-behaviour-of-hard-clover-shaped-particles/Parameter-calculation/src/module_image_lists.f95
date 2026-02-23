MODULE ImageLists

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine computes the images of a central particle.                                                                        !
!                                                                                                                                   !
! ● The subroutine takes the scaled position of a particle inside a unit box (<ScalingDistanceUnitBox>) and computes the scaled     !
!   position of the images (<ScalingDistanceImageUnitBox>).                                                                         !
!                                                                                                                                   !
! ● The images are calculated for a given number of layers (<nLayers>) in the X-, Y-, and Z-directions.                             !
! ********************************************************************************************************************************* !
SUBROUTINE ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: xLayer, yLayer, zLayer ! X-, Y-, and Z- coordinates
INTEGER( Kind= Int64 ) :: cLayer                 ! Counter (layers)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )     :: ScalingDistanceUnitBox      ! Scaled position of particles (unit box)
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( OUT ) :: ScalingDistanceImageUnitBox ! Scaled position of images (unit box)

! Position of images (unit box)
cLayer = 1
DO zLayer = - nLayers, nLayers
  DO yLayer = - nLayers, nLayers
    DO xLayer = - nLayers, nLayers
      IF( zLayer == 0 .AND. yLayer == 0 .AND. xLayer == 0 ) CYCLE ! This is not an image
      ScalingDistanceImageUnitBox(cLayer,1) = ScalingDistanceUnitBox(1) + DBLE( xLayer )
      ScalingDistanceImageUnitBox(cLayer,2) = ScalingDistanceUnitBox(2) + DBLE( yLayer )
      ScalingDistanceImageUnitBox(cLayer,3) = ScalingDistanceUnitBox(3) + DBLE( zLayer )
      cLayer = cLayer + 1
    END DO
  END DO
END DO

RETURN

END SUBROUTINE ImageConstruction

END MODULE ImageLists