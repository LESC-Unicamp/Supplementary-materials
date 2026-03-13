! ********************************************************************************************************************************* !
! This program computes order parameters, local bond parameters, structure factors, and radial distribution functions for a system  !
! of particles. The program reads a trajectory file (<FileName>) and computes the desired parameters.                               !
! ********************************************************************************************************************************* !
! Author:     Nathan Barros de Souza                                                                                                !
! Supervisor: Dr. Luís Fernando Mercier Franco and Dr. Carlos Avendaño                                                              !
! Date:       2025-03-18                                                                                                            !
! ********************************************************************************************************************************* !
PROGRAM Order_and_Bond_Parameter

! Uses eight modules: global variables, quaternions, vector operations, bond parameters, order parameters, structure assessment, variable initialization, and radial distribution functions
USE GlobalVariables
USE Quaternions
USE VectorOperations
USE BondParameters
USE OrderParameters
USE Structure
USE Initialization
USE RadialDistributionFunctions

IMPLICIT NONE

! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle             ! Counter (particle)
INTEGER( Kind= Int64 ) :: iFrame                ! Counter (frame)
INTEGER( Kind= Int64 ) :: nFrames, nEquilFrames ! Total number of frames and equilibration frames
INTEGER( Kind= Int64 ) :: iOrder                ! Counter (bond parameter)
INTEGER( Kind= Int64 ) :: iStruc, jStruc        ! Counters (structure factor)
INTEGER( Kind= Int64 ) :: ReadFrames            ! Counter (number of frames read)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: NematicOrderS2, cCubaticOrder, SmecticOrder, ChiralOrder                  ! Order parameters (nematic, cubatic, smectic, and chiral)
REAL( Kind= Real64 ) :: CubaticOrderAccumulator, NematicOrderAccumulator, SmecticOrderAccumulator ! Accumulators (order parameters)
REAL( Kind= Real64 ) :: OptimalLayerSpacing, OptimalLayerSpacingAccumulator                       ! Optimized layer spacing (smectic phase) and accumulator
REAL( Kind= Real64 ) :: ChiralAccumulator                                                         ! Accumulator (chiral order parameter)
REAL( Kind= Real64 ) :: gOrientationalOrder, gOrientationalOrderAccumulator                       ! Global orientational order parameter and accumulator

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: PhaseDirector ! Nematic phase director

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: Psi                               ! Local bond parameters
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: LocalBondParameterPsi             ! Accumulator for the local bond parameters
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: ChiralityParticle                 ! Chirality per particle
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: ChiralityParticleAccumulator      ! Chirality per particle (accumulator)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulator              ! Histogram accumulator (RDF)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulatorParallel      ! Histogram accumulator (RDF parallel to the phase director)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulatorParallelStack ! Histogram accumulator (RDF parallel to the molecular orientation)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulatorOrthogonal    ! Histogram accumulator (RDF orthogonal to the phase director)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulatorP2            ! Histogram accumulator (RDF for the second-order Legendre polynomial)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramAccumulatorP4            ! Histogram accumulator (RDF for the fourth-order Legendre polynomial)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: Histogram                         ! Histogram (RDF)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramParallel                 ! Histogram (RDF parallel to the phase director)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramParallelStack            ! Histogram (RDF parallel to the molecular orientation)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramOrthogonal               ! Histogram (RDF orthogonal to the phase director)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: HistogramP2, HistogramP4          ! Histogram (RDF for the Legendre polynomials)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: gSteinhardtOrder                  ! Steinhardt order parameter (global)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: gSteinhardtOrderAccumulator       ! Accumulator for the Steinhardt order parameter (global)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: SteinhardtOrder                   ! Steinhardt order parameter
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: SteinhardtOrderAccumulator        ! Accumulator for the Steinhardt order parameter
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: StructureFactor                   ! Structure factor
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: TotalStructureFactor              ! Structure factor (accumulated)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: gWignerOrder                      ! Wigner order parameter (global)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: gWignerOrderAccumulator           ! Accumulator for the Wigner order parameter (global)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: WignerOrder                       ! Wigner order parameter
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: WignerOrderAccumulator            ! Accumulator for the Wigner order parameter

! LOGICAL VARIABLES
LOGICAL :: FileExists = .FALSE. ! Check if a file exists

! CHARACTER STRINGS
CHARACTER( Len= 10 ) :: OrderDescriptor                ! Descriptor for the bond parameter
CHARACTER( Len= 10 ) :: SteinhardtDescriptor           ! Descriptor for the Steinhardt order parameter
CHARACTER( Len= 10 ) :: DescriptorDate, DescriptorHour ! Date and time descriptors
CHARACTER( Len= 10 ) :: FormatDate, FormatHour         ! Date and time formats
CHARACTER( Len= 70 ) :: OutputFile                     ! Output file

! Get the maximum number of threads
#ifdef _OPENMP
  nThreads = OMP_GET_MAX_THREADS(  )
#else
  nThreads = 1
#endif

! Print the number of threads
WRITE( *, "(A,I0,A)" ) "Program running with ", nThreads, " threads"
WRITE( *, "(A,I0,A)" ) " "
WRITE( *, "(A,I0,A)" ) "ATENÇÃO! Você inverteu a ordem dos quatérnios nas simulações de WXYZ para XYZW."
WRITE( *, "(A)" ) " "

! Initialize parameters
CALL Initialization_Parameters(  )

! Initialization
nFrames      = MERGE( INT( MaxCycles / nSave ), INT( DBLE( MaxCycles ) / DBLE( nSave ) ) + 1, MOD( MaxCycles, nSave ) == 0 )
nEquilFrames = MERGE( INT( ( MaxCycles - nEquilibration ) / nSave ), INT( DBLE( MaxCycles - nEquilibration ) / DBLE( nSave ) ) + &
&              1, MOD( ( MaxCycles - nEquilibration ), nSave ) == 0 )
ReadFrames = 0

! Particle dimensions
pLengthClover = pLength
pDiameterClover = 2.0D0 * pDiameter

! Allocation
ALLOCATE(  pPosition( nParticles, 3 ), pQuaternion( nParticles, 0:3 ), pOrientation( nParticles, 3 ) )
pPosition = 0.D0
pQuaternion = 0.D0
pOrientation = 0.D0
ALLOCATE( StructureFactor( -nStruc:nStruc, -nStruc:nStruc ) )
StructureFactor = 0.D0
ALLOCATE( TotalStructureFactor( -nStruc:nStruc, -nStruc:nStruc ) )
TotalStructureFactor = 0.D0
ALLOCATE( Psi( MaxNeighbors - 3 ) )
Psi = 0.D0
ALLOCATE( LocalBondParameterPsi( MaxNeighbors - 3 ) )
LocalBondParameterPsi = 0.D0
ALLOCATE( Histogram( nBins ) )
Histogram = 0.D0
ALLOCATE( HistogramParallel( nBins ) )
HistogramParallel = 0.D0
ALLOCATE( HistogramParallelStack( nBins ) )
HistogramParallelStack = 0.D0
ALLOCATE( HistogramOrthogonal( nBins ) )
HistogramOrthogonal = 0.D0
ALLOCATE( HistogramP2( nBins ) )
HistogramP2 = 0.D0
ALLOCATE( HistogramP4( nBins ) )
HistogramP4 = 0.D0
ALLOCATE( HistogramAccumulator( nBins ) )
HistogramAccumulator = 0.D0
ALLOCATE( HistogramAccumulatorParallel( nBins ) )
HistogramAccumulatorParallel = 0.D0
ALLOCATE( HistogramAccumulatorParallelStack( nBins ) )
HistogramAccumulatorParallelStack = 0.D0
ALLOCATE( HistogramAccumulatorOrthogonal( nBins ) )
HistogramAccumulatorOrthogonal = 0.D0
ALLOCATE( HistogramAccumulatorP2( nBins ) )
HistogramAccumulatorP2 = 0.D0
ALLOCATE( HistogramAccumulatorP4( nBins ) )
HistogramAccumulatorP4 = 0.D0
ALLOCATE( SteinhardtOrder( nParticles, (MaximumSteinhardt - MinimumSteinhardt) + 1 ) )
SteinhardtOrder = 0.D0
ALLOCATE( SteinhardtOrderAccumulator( nParticles, (MaximumSteinhardt - MinimumSteinhardt) + 1 ) )
SteinhardtOrderAccumulator = 0.D0
ALLOCATE( WignerOrder( nParticles, (MaximumSteinhardt - MinimumSteinhardt) + 1, 2 ) )
WignerOrder = 0.D0
ALLOCATE( WignerOrderAccumulator( nParticles, (MaximumSteinhardt - MinimumSteinhardt) + 1, 2 ) )
WignerOrderAccumulator = 0.D0
ALLOCATE( gSteinhardtOrder( (MaximumSteinhardt - MinimumSteinhardt) + 1 ) )
gSteinhardtOrder = 0.D0
ALLOCATE( gSteinhardtOrderAccumulator( (MaximumSteinhardt - MinimumSteinhardt) + 1 ) )
gSteinhardtOrderAccumulator = 0.D0
ALLOCATE( gWignerOrder( (MaximumSteinhardt - MinimumSteinhardt) + 1, 2 ) )
gWignerOrder = 0.D0
ALLOCATE( gWignerOrderAccumulator( (MaximumSteinhardt - MinimumSteinhardt) + 1, 2 ) )
gWignerOrderAccumulator = 0.D0
ALLOCATE( ChiralityParticle( nParticles ) )
ChiralityParticle = 0.D0
ALLOCATE( ChiralityParticleAccumulator( nParticles ) )
ChiralityParticleAccumulator = 0.D0

! Descriptors
WRITE( OrderDescriptor, "(I0)" ) MaxNeighbors - 3
WRITE( SteinhardtDescriptor, "(I0)" ) MaximumSteinhardt - MinimumSteinhardt + 1

! *********************************************************************************************** !
! Body-fixed position (cylinders)                                                                 !
! *********************************************************************************************** !

! First quarter
cReferencePosition(1,1) = 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(1,2) = 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(1,3) = 0.D0                               ! Å
! Second quarter
cReferencePosition(2,1) = - 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(2,2) = 0.25D0 * DSQRT( 2.D0 ) * pDiameter   ! Å
cReferencePosition(2,3) = 0.D0                                 ! Å
! Third quarter
cReferencePosition(3,1) = - 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(3,2) = - 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(3,3) = 0.D0                                 ! Å
! Fourth quarter
cReferencePosition(4,1) = 0.25D0 * DSQRT( 2.D0 ) * pDiameter   ! Å
cReferencePosition(4,2) = - 0.25D0 * DSQRT( 2.D0 ) * pDiameter ! Å
cReferencePosition(4,3) = 0.D0                                 ! Å

! Filename
WRITE( *, "(A)" ) "The filename is: "//TRIM( FileName )
INQUIRE( File= FileName, Exist= FileExists )
IF( .NOT. FileExists ) THEN
  WRITE( *, "(G0)" ) "File does not exist!"
  STOP
END IF
WRITE( *, "(A)" ) " "
WRITE( *, "(A)" ) "Starting calculation..."
WRITE( *, "(A)" ) " "

! Date and time
CALL DATE_AND_TIME( VALUES= DateTime )

! Date format (YYYY/MM/DD)
FormatDate = "(I4,2I2.2)"
WRITE( DescriptorDate, FormatDate ) DateTime(1), DateTime(2), DateTime(3)

! Time format (HH:MM:SS)
FormatHour = "(3I2.2)"
WRITE( DescriptorHour, FormatHour ) DateTime(5), DateTime(6), DateTime(7)

! Create folder to store the results
INQUIRE( File= "Calculations/", Exist= FileExists )
IF( .NOT. FileExists ) THEN
  CALL SYSTEM( "mkdir Calculations/" )
END IF

! Open file
OPEN( Unit= 10, File= FileName, Action= "READ" )

! Header of the local bond parameters
IF( ComputeLocalBondParameters ) THEN
  OutputFile = "Calculations/psi_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 15, File= OutputFile, Action= "WRITE" )
  WRITE( 15, "(G0)", Advance= "NO" ) "Frame"
  DO iOrder = 4, MaxNeighbors
    WRITE( 15, "(2G0)", Advance= "NO" ) ",ψ", iOrder
  END DO
  WRITE( 15, "(G0)", Advance= "YES" ) ",χ"
END IF

! Header of the cubatic order parameter
IF( ComputeCubaticOrder ) THEN
  OutputFile = "Calculations/cubatic_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 20, File= OutputFile, Action= "WRITE" )
  WRITE( 20, "(G0)" ) "Frame,CubaticOrder"
END IF

! Header of the nematic order parameter
IF( ComputeNematicOrder ) THEN
  OutputFile = "Calculations/nematic_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 25, File= OutputFile, Action= "WRITE" )
  WRITE( 25, "(G0)" ) "Frame,NematicOrder"
END IF

! Header of the smectic order parameter
IF( ComputeSmecticOrder ) THEN
  OutputFile = "Calculations/smectic_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 30, File= OutputFile, Action= "WRITE" )
  WRITE( 30, "(G0)" ) "Frame,SmecticOrder,OptimalLayerSpacing"
END IF

! Header of the Steinhardt order parameter
IF( ComputeSteinhardtOrder ) THEN
  OutputFile = "Calculations/steinhardt_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 35, File= OutputFile, Action= "WRITE" )
  WRITE( 35, "(G0)", Advance= "NO" ) "Particle"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 35, "(2G0)", Advance= "NO" ) ",Q", iOrder
  END DO
  WRITE( 35, "(G0)", Advance= "YES" ) ""
  OutputFile = "Calculations/steinhardt_global_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 50, File= OutputFile, Action= "WRITE" )
  WRITE( 50, "(G0)", Advance= "NO" ) "Frame"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 50, "(2G0)", Advance= "NO" ) ",Q", iOrder
  END DO
  WRITE( 50, "(G0)", Advance= "YES" ) ""
END IF

! Header of the Wigner order parameter
IF( ComputeSteinhardtOrder .AND. ComputeWignerParameter ) THEN
  OutputFile = "Calculations/wigner_real_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 40, File= OutputFile, Action= "WRITE" )
  WRITE( 40, "(G0)", Advance= "NO" ) "Particle"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 40, "(2G0)", Advance= "NO" ) ",W", iOrder
  END DO
  WRITE( 40, "(G0)", Advance= "YES" ) ""
  OutputFile = "Calculations/wigner_imag_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 45, File= OutputFile, Action= "WRITE" )
  WRITE( 45, "(G0)", Advance= "NO" ) "Particle"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 45, "(2G0)", Advance= "NO" ) ",W", iOrder
  END DO
  WRITE( 45, "(G0)", Advance= "YES" ) ""
  OutputFile = "Calculations/wigner_global_real_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 55, File= OutputFile, Action= "WRITE" )
  WRITE( 55, "(G0)", Advance= "NO" ) "Frame"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 55, "(2G0)", Advance= "NO" ) ",W", iOrder
  END DO
  WRITE( 55, "(G0)", Advance= "YES" ) ""
  OutputFile = "Calculations/wigner_global_imag_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 60, File= OutputFile, Action= "WRITE" )
  WRITE( 60, "(G0)", Advance= "NO" ) "Frame"
  DO iOrder = MinimumSteinhardt, MaximumSteinhardt
    WRITE( 60, "(2G0)", Advance= "NO" ) ",W", iOrder
  END DO
  WRITE( 60, "(G0)", Advance= "YES" ) ""
END IF

! Header of chirality per particle
IF( ComputeLocalBondParameters ) THEN
  OutputFile = "Calculations/chirality_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 65, File= OutputFile, Action= "WRITE" )
  WRITE( 65, "(G0)" ) "Particle,χ"
END IF

! Header of the global orientational order parameter
IF( ComputeGlobalOrientOrder ) THEN
  OutputFile = "Calculations/global_orientational_order_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 70, File= OutputFile, Action= "WRITE" )
  WRITE( 70, "(G0)" ) "Frame,gOrientationalOrder"
END IF

! Initialize accumulators
CubaticOrderAccumulator = 0.D0
NematicOrderAccumulator = 0.D0
SmecticOrderAccumulator = 0.D0
OptimalLayerSpacingAccumulator = 0.D0
ChiralAccumulator = 0.D0

! Loop over the frames of the trajectory file
DO iFrame = 1, nFrames
  WRITE( *, "(2A,I0,A,I0)", Advance= "NO" ) CHAR(13), "Reading frame: ", iFrame, " of ", nFrames
  READ( 10, * )
  READ( 10, * )
  READ( 10, * ) Dummy, Dummy, BoxLength(1:9) ! Box dimensions
  CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )
  DO iParticle = 1, nParticles
    READ( 10, * ) Dummy, pPosition(iParticle,1:3), pQuaternion(iParticle,0:3) ! Particle position and quaternion
    CALL VectorRotation( BodyFixedAxis, pQuaternion(iParticle,0:3), pOrientation(iParticle,1:3) ) ! Particle orientation
  END DO
  IF( iFrame <= nEquilFrames .AND. .NOT. ComputeParameterEquil ) CYCLE ! Skip equilibration frames
  IF( MOD( iFrame, nSkip + 1 ) == 0 ) THEN
    ! Compute the nematic phase director and nematic order parameter
    IF( ComputeLocalBondParameters .OR. ComputeStructureFactor .OR. ComputeRDF .OR. ComputeNematicOrder &
    &   .OR. ComputeSmecticOrder .OR. ComputeGlobalOrientOrder ) THEN
      CALL NematicOrderParameter( pOrientation, NematicOrderS2, PhaseDirector )
      IF( ComputeNematicOrder ) THEN
        NematicOrderAccumulator = NematicOrderAccumulator + NematicOrderS2
        WRITE( 25, "(I0,', ',G0.10)" ) iFrame, NematicOrderS2
        FLUSH( 25 )
      END IF
    END IF
    ! Compute the local bond parameters
    IF( ComputeLocalBondParameters ) THEN
      CALL LocalBondParameters( PhaseDirector, Psi, ChiralOrder, ChiralityParticle )
      LocalBondParameterPsi = LocalBondParameterPsi + Psi
      ChiralAccumulator = ChiralAccumulator + ChiralOrder
      ChiralityParticleAccumulator = ChiralityParticleAccumulator + ChiralityParticle
      WRITE( 15, "(I0,"//TRIM( OrderDescriptor )//"(',',G0.10),',',G0.10)" ) iFrame, Psi, ChiralOrder
      FLUSH( 15 )
    END IF
    ! Compute the structure factor
    IF( ComputeStructureFactor .OR. ComputeGlobalOrientOrder ) CALL Structure_Factor( iFrame, nFrames, nStruc, PhaseDirector, &
    &                                                                                 StructureFactor, gOrientationalOrder )
    IF( ComputeStructureFactor ) TotalStructureFactor = TotalStructureFactor + StructureFactor
    IF( ComputeGlobalOrientOrder ) THEN
      gOrientationalOrderAccumulator = gOrientationalOrderAccumulator + gOrientationalOrder
      WRITE( 70, "(I0,', ',G0.10)" ) iFrame, gOrientationalOrder
      FLUSH( 70 )
    END IF
    ! Compute the radial distribution functions
    IF( ComputeRDF ) THEN
      CALL RDFunction( HistogramAccumulator ) ! RDF
      Histogram = Histogram + HistogramAccumulator
      CALL RDFunction_Parallel( PhaseDirector, HistogramAccumulatorParallel ) ! RDF parallel to the phase director
      HistogramParallel = HistogramParallel + HistogramAccumulatorParallel
      CALL RDFunction_Parallel_Stack( pOrientation, HistogramAccumulatorParallelStack ) ! RDF parallel to the molecular orientation
      HistogramParallelStack = HistogramParallelStack + HistogramAccumulatorParallelStack
      CALL RDFunction_Orthogonal( PhaseDirector, HistogramAccumulatorOrthogonal ) ! RDF orthogonal to the phase director
      HistogramOrthogonal = HistogramOrthogonal + HistogramAccumulatorOrthogonal
      CALL RDFunction_LegendrePolynomials( pOrientation, HistogramAccumulatorP2, HistogramAccumulatorP4 ) ! RDF for the Legendre polynomials
      HistogramP2 = HistogramP2 + HistogramAccumulatorP2
      HistogramP4 = HistogramP4 + HistogramAccumulatorP4
    END IF
    ! Compute the cubatic order parameter
    IF( ComputeCubaticOrder ) THEN
      CALL CubaticOrderParameter( iFrame, nFrames, pQuaternion, cCubaticOrder )
      CubaticOrderAccumulator = CubaticOrderAccumulator + cCubaticOrder
      WRITE( 20, "(I0,', ',G0.10)" ) iFrame, cCubaticOrder
      FLUSH( 20 )
    END IF
    ! Compute the smectic order parameter
    IF( ComputeSmecticOrder ) THEN
      CALL SmecticOrderParameter( iFrame, nFrames, pPosition, PhaseDirector, SmecticOrder, OptimalLayerSpacing )
      SmecticOrderAccumulator = SmecticOrderAccumulator + SmecticOrder
      OptimalLayerSpacingAccumulator = OptimalLayerSpacingAccumulator + OptimalLayerSpacing
      WRITE( 30, "(I0,', ',G0.10,', ',G0.10)" ) iFrame, SmecticOrder, OptimalLayerSpacing
      FLUSH( 30 )
    END IF
    ! Compute the rotational invariants (Steinhardt and Wigner order parameters)
    IF( ComputeSteinhardtOrder ) THEN
      CALL Steinhardt_Parameters( iFrame, nFrames, SteinhardtOrder, WignerOrder, gSteinhardtOrder, gWignerOrder )
      DO iOrder = 1, MaximumSteinhardt - MinimumSteinhardt + 1
        SteinhardtOrderAccumulator(:,iOrder) = SteinhardtOrderAccumulator(:,iOrder) + SteinhardtOrder(:,iOrder)
        gSteinhardtOrderAccumulator(iOrder) = gSteinhardtOrderAccumulator(iOrder) + gSteinhardtOrder(iOrder)
        IF( ComputeWignerParameter ) THEN
          WignerOrderAccumulator(:,iOrder,1) = WignerOrderAccumulator(:,iOrder,1) + WignerOrder(:,iOrder,1)
          WignerOrderAccumulator(:,iOrder,2) = WignerOrderAccumulator(:,iOrder,2) + WignerOrder(:,iOrder,2)
          gWignerOrderAccumulator(iOrder,1) = gWignerOrderAccumulator(iOrder,1) + gWignerOrder(iOrder,1)
          gWignerOrderAccumulator(iOrder,2) = gWignerOrderAccumulator(iOrder,2) + gWignerOrder(iOrder,2)
        END IF
      END DO
      WRITE( 50, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iFrame, gSteinhardtOrder
      WRITE( 55, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iFrame, gWignerOrder(:,1)
      WRITE( 60, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iFrame, gWignerOrder(:,2)
    END IF
    ! Reset
    pQuaternion = 0.D0
    BoxLength   = 0.D0
    pPosition   = 0.D0
    ! Update the number of frames read
    ReadFrames = ReadFrames + 1
  END IF
END DO

! Close trajectory file
CLOSE( 10 )

! Status
WRITE( *, "(A)" ) " "
WRITE( *, "(A)" ) " "

! Average the results
IF( ComputeLocalBondParameters ) LocalBondParameterPsi = LocalBondParameterPsi / DBLE( ReadFrames )
IF( ComputeLocalBondParameters ) ChiralAccumulator = ChiralAccumulator / DBLE( ReadFrames )
IF( ComputeLocalBondParameters ) ChiralityParticleAccumulator = ChiralityParticleAccumulator / DBLE( ReadFrames )
IF( ComputeStructureFactor ) TotalStructureFactor = TotalStructureFactor / DBLE( ReadFrames )
IF( ComputeRDF ) THEN
  Histogram = Histogram / DBLE( ReadFrames )
  HistogramParallel = HistogramParallel / DBLE( ReadFrames )
  HistogramParallelStack = HistogramParallelStack / DBLE( ReadFrames )
  HistogramOrthogonal = HistogramOrthogonal / DBLE( ReadFrames )
  HistogramP2 = HistogramP2 / DBLE( ReadFrames )
  HistogramP4 = HistogramP4 / DBLE( ReadFrames )
END IF
IF( ComputeCubaticOrder ) THEN
  CubaticOrderAccumulator = CubaticOrderAccumulator / DBLE( ReadFrames )
  WRITE( 20, "(A)" ) " "
  WRITE( *, "(A,G0.10)" ) "Cubatic order parameter: ", CubaticOrderAccumulator
  WRITE( 20, "(A,G0.10)" ) "Average cubatic order parameter: ", CubaticOrderAccumulator
  FLUSH( 20 )
  WRITE( *, "(A)" ) " "
END IF
IF( ComputeGlobalOrientOrder ) THEN
  gOrientationalOrderAccumulator = gOrientationalOrderAccumulator / DBLE( ReadFrames )
  WRITE( 70, "(A)" ) " "
  WRITE( *, "(A,G0.10)" ) "Global orientational order parameter: ", gOrientationalOrderAccumulator
  WRITE( 70, "(A,G0.10)" ) "Average global orientational order parameter: ", gOrientationalOrderAccumulator
  FLUSH( 70 )
  WRITE( *, "(A)" ) " "
END IF
IF( ComputeNematicOrder ) THEN
  NematicOrderAccumulator = NematicOrderAccumulator / DBLE( ReadFrames )
  WRITE( 25, "(A)" ) " "
  WRITE( *, "(A,G0.10)" ) "Nematic order parameter: ", NematicOrderAccumulator
  WRITE( 25, "(A,G0.10)" ) "Average nematic order parameter: ", NematicOrderAccumulator
  FLUSH( 25 )
  WRITE( *, "(A)" ) " "
END IF
IF( ComputeSmecticOrder ) THEN
  SmecticOrderAccumulator = SmecticOrderAccumulator / DBLE( ReadFrames )
  OptimalLayerSpacingAccumulator = OptimalLayerSpacingAccumulator / DBLE( ReadFrames )
  WRITE( 30, "(A)" ) " "
  WRITE( *, "(A,G0.10)" ) "Smectic order parameter: ", SmecticOrderAccumulator
  WRITE( *, "(A,G0.10)" ) "Optimal layer spacing: ", OptimalLayerSpacingAccumulator
  WRITE( 30, "(A,G0.10)" ) "Average smectic order parameter: ", SmecticOrderAccumulator
  WRITE( 30, "(A,G0.10)" ) "Average optimal layer spacing: ", OptimalLayerSpacingAccumulator
  FLUSH( 30 )
  WRITE( *, "(A)" ) " "
END IF
IF( ComputeSteinhardtOrder ) THEN
  SteinhardtOrderAccumulator = SteinhardtOrderAccumulator / DBLE( ReadFrames )
  gSteinhardtOrderAccumulator = gSteinhardtOrderAccumulator / DBLE( ReadFrames )
  IF( ComputeWignerParameter ) WignerOrderAccumulator = WignerOrderAccumulator / DBLE( ReadFrames )
  IF( ComputeWignerParameter ) gWignerOrderAccumulator = gWignerOrderAccumulator / DBLE( ReadFrames )
  DO iParticle = 1, nParticles
    WRITE( 35, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iParticle, SteinhardtOrderAccumulator(iParticle,:)
    IF( ComputeWignerParameter ) WRITE( 40, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iParticle, &
    &                            WignerOrderAccumulator(iParticle,:,1)
    IF( ComputeWignerParameter ) WRITE( 45, "(I0,"//TRIM( SteinhardtDescriptor )//"(',',G0.10))" ) iParticle, &
    &                            WignerOrderAccumulator(iParticle,:,2)
  END DO
  FLUSH( 35 )
  IF( ComputeWignerParameter ) FLUSH( 40 )
  IF( ComputeWignerParameter ) FLUSH( 45 )
END IF

! Chirality per particle
IF( ComputeLocalBondParameters ) THEN
  DO iParticle = 1, nParticles
    WRITE( 65, "(I0,',',G0.10)" ) iParticle, ChiralityParticleAccumulator(iParticle)
  END DO
  FLUSH( 65 )
END IF

! Cap the structure factor
IF( ApplyCapping .AND. ComputeLocalBondParameters ) THEN
  WHERE( TotalStructureFactor <= CappingFactor ) TotalStructureFactor = 0.D0
  WHERE( TotalStructureFactor > CappingFactor  ) TotalStructureFactor = 1.D0
END IF

IF( ComputeLocalBondParameters ) THEN
  WRITE( *, "(A,F8.6,A)" ) "Bond Parameter ψ4: ", LocalBondParameterPsi(1)
  WRITE( *, "(A,F8.6,A)" ) "Bond Parameter ψ5: ", LocalBondParameterPsi(2)
  WRITE( *, "(A,F8.6,A)" ) "Bond Parameter ψ6: ", LocalBondParameterPsi(3)
  WRITE( *, "(A,F8.6,A)" ) "Chiral order parameter: ", ChiralAccumulator
  WRITE( *, "(A)" ) " "
  WRITE( 15, "(A)" ) " "
  WRITE( 15, "(A,G0.10)" ) "Average bond parameter ψ4: ", LocalBondParameterPsi(1)
  WRITE( 15, "(A,G0.10)" ) "Average bond parameter ψ5: ", LocalBondParameterPsi(2)
  WRITE( 15, "(A,G0.10)" ) "Average bond parameter ψ6: ", LocalBondParameterPsi(3)
  WRITE( 15, "(A,G0.10)" ) "Average chiral order parameter: ", ChiralAccumulator
END IF

IF( ComputeSteinhardtOrder ) THEN
  DO iOrder = 1, MaximumSteinhardt - MinimumSteinhardt + 1
    WRITE( *, "(A,G0.10)" ) "Global Steinhardt order parameter: ", gSteinhardtOrderAccumulator(iOrder)
    WRITE( *, "(A,G0.10)" ) "Global Wigner order parameter (real): ", gWignerOrderAccumulator(iOrder,1)
    WRITE( *, "(A,G0.10)" ) "Global Wigner order parameter (imaginary): ", gWignerOrderAccumulator(iOrder,2)
    WRITE( *, "(A)" ) " "
  END DO
  WRITE( 50, "(A)" ) " "
  WRITE( 55, "(A)" ) " "
  WRITE( 60, "(A)" ) " "
  DO iOrder = 1, MaximumSteinhardt - MinimumSteinhardt + 1
    WRITE( 50, "(A,G0.10)" ) "Average global Steinhardt order parameter: ", gSteinhardtOrderAccumulator(iOrder)
    WRITE( 55, "(A,G0.10)" ) "Average global Wigner order parameter (real): ", gWignerOrderAccumulator(iOrder,1)
    WRITE( 60, "(A,G0.10)" ) "Average global Wigner order parameter (imaginary): ", gWignerOrderAccumulator(iOrder,2)
  END DO
END IF

! Close files
IF( ComputeLocalBondParameters ) CLOSE( 15 )
IF( ComputeCubaticOrder ) CLOSE( 20 )
IF( ComputeNematicOrder ) CLOSE( 25 )
IF( ComputeSmecticOrder ) CLOSE( 30 )
IF( ComputeSteinhardtOrder ) CLOSE( 35 )
IF( ComputeWignerParameter ) CLOSE( 40 )
IF( ComputeWignerParameter ) CLOSE( 45 )
IF( ComputeGlobalOrientOrder ) CLOSE( 70 )

! Write out histograms in terms of the radial distance
IF( ComputeRDF ) THEN
  ! RDF
  OutputFile = "Calculations/gr_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "RadialDistance,RadialDistributionFunction"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, Histogram(iParticle)
  END DO
  CLOSE( 10 )
  ! RDF parallel to the phase director
  OutputFile = "Calculations/gr_parallel_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "ParallelRadialDistance,RadialDistributionFunction"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, HistogramParallel(iParticle)
  END DO
  CLOSE( 10 )
  ! RDF parallel to the molecular orientation
  OutputFile = "Calculations/gr_parallel_stacks_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "ParallelRadialDistance,RadialDistributionFunction"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, HistogramParallelStack(iParticle)
  END DO
  CLOSE( 10 )
  ! RDF orthogonal to the phase director
  OutputFile = "Calculations/gr_orthogonal_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "OrthogonalRadialDistance,RadialDistributionFunction"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, HistogramOrthogonal(iParticle)
  END DO
  CLOSE( 10 )
  ! RDF for the second-order Legendre polynomial
  OutputFile = "Calculations/gr_p2_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "RadialDistance,SecondLegendrePolynomial"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, HistogramP2(iParticle)
  END DO
  CLOSE( 10 )
  ! RDF for the fourth-order Legendre polynomial
  OutputFile = "Calculations/gr_p4_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "RadialDistance,FourthLegendrePolynomial"
  DO iParticle = 1, nBins
    WRITE( 10, "(G0.10,',',G0.10)" ) DBLE( iParticle - 1 ) * bWidth, HistogramP4(iParticle)
  END DO
  CLOSE( 10 )
END IF

! Write out the structure factor
IF( ComputeStructureFactor ) THEN
  OutputFile = "Calculations/structure_factor_"//TRIM( DescriptorDate )//"_"//TRIM( DescriptorHour )//".dat"
  OPEN( Unit= 10, File= OutputFile )
  WRITE( 10, "(A)" ) "x,y,S(k)"
  DO iStruc = -nStruc, nStruc
    DO jStruc = -nStruc, nStruc
      IF( ApplyCapping ) WRITE( 10, "(2(I0,','),I0)" ) iStruc, jStruc, INT( ANINT( TotalStructureFactor(iStruc,jStruc) ) )
      IF( .NOT. ApplyCapping ) WRITE( 10, "(2(I0,','),G0.10)" ) iStruc, jStruc, TotalStructureFactor(iStruc,jStruc)
    END DO
  END DO
  CLOSE( 10 )
END IF

! Status
WRITE( *, "(A)" ) "Calculation completed!"

END PROGRAM Order_and_Bond_Parameter