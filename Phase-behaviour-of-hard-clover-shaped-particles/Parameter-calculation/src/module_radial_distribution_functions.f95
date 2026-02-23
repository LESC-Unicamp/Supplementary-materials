MODULE RadialDistributionFunctions

! Uses three modules: global variables, vector operations, and quaternions
USE GlobalVariables
USE VectorOperations
USE Quaternions

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the radial distribution function (<Histogram>) of the system.                                          !
! ********************************************************************************************************************************* !
SUBROUTINE RDFunction( Histogram )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Particle counters
INTEGER( Kind= Int64 ) :: bBin                 ! Bin counter

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Distance        ! Distance
REAL( Kind= Real64 ) :: Density         ! Number density
REAL( Kind= Real64 ) :: bVolume         ! Volume of the bin
REAL( Kind= Real64 ) :: bPreviousVolume ! Volume of the bin (previous bin)
REAL( Kind= Real64 ) :: bNextVolume     ! Volume of the bin (next bin)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                :: VectorDistance         ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                :: ScalingDistanceUnitBox ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT ) :: Histogram              ! Histogram of the radial distribution function

! Initialization
Histogram = 0.D0

! Calculate the number density
Density = nParticles / BoxVolume

! Loop over all particle pairs
DO iParticle = 1, nParticles - 1
  DO jParticle = iParticle + 1, nParticles
    ! Calculate distance between particles i and j
    VectorDistance(1:3) = pPosition(iParticle,1:3) - pPosition(jParticle,1:3)
    ! Apply minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Calculate the distance
    Distance = DOT_PRODUCT( VectorDistance, VectorDistance )
    Distance = DSQRT( Distance )
    ! Assign the distance to the corresponding bin
    bBin = FLOOR( Distance / bWidth ) + 1
    IF (bBin <= nBins) THEN
      Histogram(bBin) = Histogram(bBin) + 2.D0 ! Each particle pair is counted twice: (i,j) and (j,i)
    END IF
  END DO
END DO

! Normalize the histogram
DO bBin = 1, nBins
  bPreviousVolume = ( 4.D0 / 3.D0 ) * cPi * ( ( bBin - 1 ) * bWidth ) * ( ( bBin - 1 ) * bWidth ) * ( ( bBin - 1 ) * bWidth ) ! Volume of the sphere (previous bin)
  bNextVolume = ( 4.D0 / 3.D0 ) * cPi * ( bBin * bWidth ) * ( bBin * bWidth ) * ( bBin * bWidth ) ! Volume of the sphere (current bin)
  bVolume = bNextVolume - bPreviousVolume ! Volume of the current bin
  Histogram(bBin) = Histogram(bBin) / ( Density * nParticles * bVolume ) ! Normalization
END DO

RETURN

END SUBROUTINE RDFunction

! ********************************************************************************************************************************* !
! This subroutine calculates the parallel radial distribution function (<Histogram>) of the system.                                 !
!                                                                                                                                   !
! ● The histogram considers the vector distance between particles projected onto the phase director (<PhaseDirector>).              !
!                                                                                                                                   !
! ● The condition for the particles to be considered in the same column (bin) is that the orthogonal distance is less than half the !
!   diameter of the particles.                                                                                                      !
! ********************************************************************************************************************************* !
! For more information, please read: Discotic columnar liquid crystal studied in the bulk and nanoconfined states by molecular      !
!                                    dynamics simulation, J. Chem. Phys. 141, 134902 (2014), DOI: 10.1063/1.4896052.                !
! ********************************************************************************************************************************* !
SUBROUTINE RDFunction_Parallel( PhaseDirector, Histogram )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Particle counters
INTEGER( Kind= Int64 ) :: bBin                 ! Bin counter

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ParallelDistance   ! Distance (parallel)
REAL( Kind= Real64 ) :: OrthogonalDistance ! Distance (orthogonal)
REAL( Kind= Real64 ) :: Density            ! Number density
REAL( Kind= Real64 ) :: bVolume            ! Volume of the bin
REAL( Kind= Real64 ) :: bPreviousVolume    ! Volume of the bin (previous)
REAL( Kind= Real64 ) :: bNextVolume        ! Volume of the bin (next)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                :: VectorDistance           ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                :: ParallelVectorDistance   ! Vector distance projected onto the phase director (parallel)
REAL( Kind= Real64 ), DIMENSION( 3 )                :: OrthogonalVectorDistance ! Vector distance orthogonal to the phase director
REAL( Kind= Real64 ), DIMENSION( 3 )                :: ScalingDistanceUnitBox   ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )  :: PhaseDirector            ! Phase director (nematic director)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT ) :: Histogram                ! Histogram of the radial distribution function

! Initialization
Histogram = 0.D0

! Calculate the number density
Density = nParticles / BoxVolume

! Loop over all particle pairs
DO iParticle = 1, nParticles - 1
  DO jParticle = iParticle + 1, nParticles
    ! Calculate distance between particles i and j
    VectorDistance(1:3) = pPosition(iParticle,1:3) - pPosition(jParticle,1:3)
    ! Apply minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Project the distance onto the phase director (parallel)
    ParallelVectorDistance = DOT_PRODUCT( VectorDistance, PhaseDirector ) * PhaseDirector
    ! Perpendicular vector
    OrthogonalVectorDistance = VectorDistance - ParallelVectorDistance
    ! Calculate the distance (parallel)
    ParallelDistance = DOT_PRODUCT( ParallelVectorDistance, ParallelVectorDistance )
    ParallelDistance = DSQRT( ParallelDistance )
    ! Calculate the distance (orthogonal)
    OrthogonalDistance = DOT_PRODUCT( OrthogonalVectorDistance, OrthogonalVectorDistance )
    OrthogonalDistance = DSQRT( OrthogonalDistance )
    ! Condition (particles lie in the same column)
    IF( OrthogonalDistance < 0.5D0 * pDiameterClover * ParallelRange ) THEN
      ! Assign the distance to the corresponding bin
      bBin = FLOOR( ParallelDistance / bWidth ) + 1
      IF (bBin <= nBins) THEN
        Histogram(bBin) = Histogram(bBin) + 2.D0 ! Each particle pair is counted twice: (i,j) and (j,i)
      END IF
    END IF
  END DO
END DO

! Normalize the histogram
DO bBin = 1, nBins
  bPreviousVolume = 0.5D0 * cPi * pDiameterClover * pDiameterClover * ( ( bBin - 1 ) * bWidth ) * ParallelRange * ParallelRange ! Volume of the cylindrical region (previous bin): Area of the clover times twice the width of the bin (consider the region below and above the bin)
  bNextVolume = 0.5D0 * cPi * pDiameterClover * pDiameterClover * ( bBin * bWidth ) * ParallelRange * ParallelRange ! Volume of the cylindrical region (current bin): Area of the clover times twice the width of the bin (consider the region below and above the bin)
  bVolume = bNextVolume - bPreviousVolume ! Volume of the current bin
  Histogram(bBin) = Histogram(bBin) / ( Density * nParticles * bVolume ) ! Normalization
END DO

RETURN

END SUBROUTINE RDFunction_Parallel

! ********************************************************************************************************************************* !
!   This subroutine calculates the orthogonal radial distribution function (<Histogram>) of the system.                             !
!                                                                                                                                   !
! ● The histogram considers the vector distance between particles projected orthogonally to the phase director (<PhaseDirector>).   !
!                                                                                                                                   !
! ● The condition for the particles to be considered in the same layer (bin) is that the parallel distance is less than half the    !
!   length of the particles.                                                                                                        !
! ********************************************************************************************************************************* !
! For more information, please read: Discotic columnar liquid crystal studied in the bulk and nanoconfined states by molecular      !
!                                    dynamics simulation, J. Chem. Phys. 141, 134902 (2014), DOI: 10.1063/1.4896052.                !
! ********************************************************************************************************************************* !
SUBROUTINE RDFunction_Orthogonal( PhaseDirector, Histogram )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Particle counters
INTEGER( Kind= Int64 ) :: bBin                 ! Bin counter

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ParallelDistance   ! Distance (parallel)
REAL( Kind= Real64 ) :: OrthogonalDistance ! Distance (orthogonal)
REAL( Kind= Real64 ) :: Density            ! Number density
REAL( Kind= Real64 ) :: bVolume            ! Volume of the bin
REAL( Kind= Real64 ) :: bPreviousVolume    ! Volume of the bin (previous)
REAL( Kind= Real64 ) :: bNextVolume        ! Volume of the bin (next)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                :: VectorDistance           ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                :: ParallelVectorDistance   ! Vector distance projected onto the phase director (parallel)
REAL( Kind= Real64 ), DIMENSION( 3 )                :: OrthogonalVectorDistance ! Vector distance orthogonal to the phase director
REAL( Kind= Real64 ), DIMENSION( 3 )                :: ScalingDistanceUnitBox   ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )  :: PhaseDirector            ! Phase director (nematic director)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT ) :: Histogram                ! Histogram of the radial distribution function

! Initialization
Histogram = 0.D0

! Calculate the number density
Density = nParticles / BoxVolume

! Loop over all particle pairs
DO iParticle = 1, nParticles - 1
  DO jParticle = iParticle + 1, nParticles
    ! Calculate distance between particles i and j
    VectorDistance(1:3) = pPosition(iParticle,1:3) - pPosition(jParticle,1:3)
    ! Apply minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Project the distance onto the phase director (parallel)
    ParallelVectorDistance = DOT_PRODUCT( VectorDistance, PhaseDirector ) * PhaseDirector
    ! Perpendicular vector
    OrthogonalVectorDistance = VectorDistance - ParallelVectorDistance
    ! Calculate the distance (parallel)
    ParallelDistance = DOT_PRODUCT( ParallelVectorDistance, ParallelVectorDistance )
    ParallelDistance = DSQRT( ParallelDistance )
    ! Calculate the distance (orthogonal)
    OrthogonalDistance = DOT_PRODUCT( OrthogonalVectorDistance, OrthogonalVectorDistance )
    OrthogonalDistance = DSQRT( OrthogonalDistance )
    ! Condition (particles lie in the same layer)
    IF( ParallelDistance < 0.5D0 * pLengthClover * OrthogonalRange ) THEN
      ! Assign the distance to the corresponding bin
      bBin = FLOOR( OrthogonalDistance / bWidth ) + 1
      IF (bBin <= nBins) THEN
        Histogram(bBin) = Histogram(bBin) + 2.D0 ! Each particle pair is counted twice: (i,j) and (j,i)
      END IF
    END IF
  END DO
END DO

! Normalize the histogram
DO bBin = 1, nBins
  bPreviousVolume = cPi * pLengthClover * ( ( bBin - 1 ) * bWidth ) * ( ( bBin - 1 ) * bWidth ) * OrthogonalRange ! Volume of the cylindrical region (previous bin): πr²L
  bNextVolume = cPi * pLengthClover * ( bBin * bWidth ) * ( bBin * bWidth ) * OrthogonalRange ! Volume of the cylindrical region (current bin): πr²L
  bVolume = bNextVolume - bPreviousVolume ! Volume of the current bin
  Histogram(bBin) = Histogram(bBin) / ( Density * nParticles * bVolume ) ! Normalization
END DO

RETURN

END SUBROUTINE RDFunction_Orthogonal

! ********************************************************************************************************************************* !
! This subroutine calculates the parallel radial distribution function (<Histogram>) of a stack of particles.                       !
!                                                                                                                                   !
! ● The histogram considers the vector distance between particles projected onto the orientation (<pCurrentOrientation>) of the     !
!   particles.                                                                                                                      !
! ********************************************************************************************************************************* !
SUBROUTINE RDFunction_Parallel_Stack( pCurrentOrientation, Histogram )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Particle counters
INTEGER( Kind= Int64 ) :: bBin                 ! Bin counter

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ParallelDistance   ! Distance (parallel)
REAL( Kind= Real64 ) :: OrthogonalDistance ! Distance (orthogonal)
REAL( Kind= Real64 ) :: Density            ! Number density
REAL( Kind= Real64 ) :: bVolume            ! Volume of the bin
REAL( Kind= Real64 ) :: bPreviousVolume    ! Volume of the bin (previous)
REAL( Kind= Real64 ) :: bNextVolume        ! Volume of the bin (next)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: VectorDistance           ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: ParallelVectorDistance   ! Vector distance projected onto the phase director (parallel)
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: OrthogonalVectorDistance ! Vector distance orthogonal to the phase director
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: ScalingDistanceUnitBox   ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT )   :: Histogram                ! Histogram of the radial distribution function
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN ) :: pCurrentOrientation      ! Orientation of the particles

! Initialization
Histogram = 0.D0

! Calculate the number density
Density = nParticles / BoxVolume

! Loop over all particle pairs
DO iParticle = 1, nParticles
  DO jParticle = 1, nParticles
    ! Cycle if the particles are the same
    IF( iParticle == jParticle ) CYCLE
    ! Calculate distance between particles i and j
    VectorDistance(1:3) = pPosition(iParticle,1:3) - pPosition(jParticle,1:3)
    ! Apply minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Project the distance onto the orientation of particle i (parallel)
    ParallelVectorDistance = DOT_PRODUCT( VectorDistance, pCurrentOrientation(iParticle,1:3) ) * pCurrentOrientation(iParticle,1:3)
    ! Perpendicular vector
    OrthogonalVectorDistance = VectorDistance - ParallelVectorDistance
    ! Calculate the distance (parallel)
    ParallelDistance = DOT_PRODUCT( ParallelVectorDistance, ParallelVectorDistance )
    ParallelDistance = DSQRT( ParallelDistance )
    ! Calculate the distance (orthogonal)
    OrthogonalDistance = DOT_PRODUCT( OrthogonalVectorDistance, OrthogonalVectorDistance )
    OrthogonalDistance = DSQRT( OrthogonalDistance )
    ! Condition (particles lie in the same column)
    IF( OrthogonalDistance < 0.5D0 * pDiameterClover ) THEN
      ! Assign the distance to the corresponding bin
      bBin = FLOOR( ParallelDistance / bWidth ) + 1
      IF (bBin <= nBins) THEN
        Histogram(bBin) = Histogram(bBin) + 1.D0 ! Each particle pair is counted once
      END IF
    END IF
  END DO
END DO

! Normalize the histogram
DO bBin = 1, nBins
  bPreviousVolume = 0.5D0 * cPi * pDiameterClover * pDiameterClover * ( ( bBin - 1 ) * bWidth ) ! Volume of the cylindrical region (previous bin): Area of the clover times twice the width of the bin (consider the region below and above the bin)
  bNextVolume = 0.5D0 * cPi * pDiameterClover * pDiameterClover * ( bBin * bWidth ) ! Volume of the cylindrical region (current bin): Area of the clover times twice the width of the bin (consider the region below and above the bin)
  bVolume = bNextVolume - bPreviousVolume ! Volume of the current bin
  Histogram(bBin) = Histogram(bBin) / ( Density * nParticles * bVolume ) ! Normalization
END DO

RETURN

END SUBROUTINE RDFunction_Parallel_Stack

! ********************************************************************************************************************************* !
! This subroutine calculates the radial distribution function of the second (<HistogramS2>) and fourth (<HistogramS4>) Legendre     !
! polynomials, which depend on the orientation (<pCurrentOrientation>) of the particles.                                            !
! ********************************************************************************************************************************* !
SUBROUTINE RDFunction_LegendrePolynomials( pCurrentOrientation, HistogramS2, HistogramS4 )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Particle counters
INTEGER( Kind= Int64 ) :: bBin                 ! Bin counter

! REAL VARIABLES (SCALAR,CONSTANT)
REAL( Kind= Real64 ), PARAMETER :: Epsilon = 1.D-14 ! Floatting point tolerance

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Distance ! Distance

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: VectorDistance           ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: ScalingDistanceUnitBox   ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT )   :: HistogramS2, HistogramS4 ! Histogram of the radial distribution function of the second and fourth Legendre polynomials
REAL( Kind= Real64 ), DIMENSION( nBins )              :: Histogram                ! Histogram of the radial distribution function
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN ) :: pCurrentOrientation      ! Orientation of the particles

! Initialization
HistogramS2 = 0.D0
HistogramS4 = 0.D0
Histogram   = 0.D0

! Loop over all particle pairs
DO iParticle = 1, nParticles - 1
  DO jParticle = iParticle + 1, nParticles
    ! Cycle if the particles are the same
    IF( iParticle == jParticle ) CYCLE
    ! Calculate distance between particles i and j
    VectorDistance(1:3) = pPosition(iParticle,1:3) - pPosition(jParticle,1:3)
    ! Apply minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Calculate the distance
    Distance = DOT_PRODUCT( VectorDistance, VectorDistance )
    Distance = DSQRT( Distance )
    ! Assign the distance to the corresponding bin
    bBin = FLOOR( Distance / bWidth ) + 1
    IF (bBin <= nBins) THEN
      Histogram(bBin) = Histogram(bBin) + 2.D0 ! Count every particle in the bin
      HistogramS2(bBin) = HistogramS2(bBin) + 3.D0 * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) - 1.D0 ! Second Legendre polynomial (x2)
      HistogramS4(bBin) = HistogramS4(bBin) + ( 35.D0 * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) - 30.D0 * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) * DOT_PRODUCT( pCurrentOrientation(iParticle,1:3), &
      &                   pCurrentOrientation(jParticle,1:3) ) + 3.D0 ) * 0.25D0 ! Fourth Legendre polynomial (x2)
    END IF
  END DO
END DO

! Normalize the histogram
DO bBin = 1, nBins
  IF( DABS( Histogram(bBin) ) <= Epsilon ) THEN ! Avoid division by zero
    HistogramS2(bBin) = 0.D0
    HistogramS4(bBin) = 0.D0
  ELSE
    HistogramS2(bBin) = HistogramS2(bBin) / ( Histogram(bBin) + Epsilon ) ! Normalization of the second Legendre polynomial with respect to the number of particles in the bin
    HistogramS4(bBin) = HistogramS4(bBin) / ( Histogram(bBin) + Epsilon ) ! Normalization of the fourth Legendre polynomial with respect to the number of particles in the bin
  END IF
END DO

RETURN

END SUBROUTINE RDFunction_LegendrePolynomials

END MODULE RadialDistributionFunctions