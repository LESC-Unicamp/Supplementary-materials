MODULE Neighbors

! Uses 1 module: global variables
USE GlobalVariables

CONTAINS

! ********************************************************************************************************************************* !
! Subroutine to calculate the Voronoi neighbors of a given particle in 2D (plane).                                                  !
! ********************************************************************************************************************************* !
! See Allen and Tildesley, Computer Simulation of Liquids, 1987, for more information or visit                                      !
! <https://server.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/>.                                                      !
! ********************************************************************************************************************************* !
SUBROUTINE VoronoiConstruction2D( nCandidates, nVertices, nEdges, cRelativePosition, cSquaredRelativeDistance, cEdges, &
&                                 vRelativePosition, vIndex )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )                :: iCand, jCand  ! Counters (candidates)
INTEGER( Kind= Int64 )                :: lCand         ! Counter (test candidates)
INTEGER( Kind= Int64 )                :: ValidVertices ! Number of valid vertices
INTEGER( Kind= Int64 )                :: iVertex       ! Counters (vertices)
INTEGER( Kind= Int64 ), INTENT( IN )  :: nCandidates   ! Number of candidates
INTEGER( Kind= Int64 ), INTENT( OUT ) :: nVertices     ! Number of vertices
INTEGER( Kind= Int64 ), INTENT( OUT ) :: nEdges        ! Number of edges

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), INTENT( OUT ), DIMENSION( : )    :: cEdges ! Edges of the candidate
INTEGER( Kind= Int64 ), INTENT( OUT ), DIMENSION( :, : ) :: vIndex ! Index of the vertices

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Determinant        ! Determinant (check if the planes intersect)
REAL( Kind= Real64 ) :: DeterminantInverse ! Inverse of determinant

! REAL VARIABLES (SCALAR,PARAMETER)
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-10 ! Tolerance

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 2 )                   :: vCoordinates             ! Vertex coordinates
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )     :: cSquaredRelativeDistance ! Squared relative distance
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN )  :: cRelativePosition        ! Relative position of the candidates
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( OUT ) :: vRelativePosition        ! Relative position of the vertices

! LOGICAL VARIABLES
LOGICAL :: VertexValidity = .FALSE.

! Initialization
ValidVertices = 0

! Examination of each possible vertex defined by the intersection of 2 edges
DO iCand = 1, nCandidates - 1

  DO jCand = iCand + 1, nCandidates

    ! Calculate the determinant and check if the edges intersect
    Determinant = cRelativePosition(iCand,1) * cRelativePosition(jCand,2) - cRelativePosition(jCand,1) * cRelativePosition(iCand,2)

    ! Condition for edge intersection (non-zero determinant means unique solution or a unique point)
    IF( DABS( Determinant ) > Tolerance ) THEN

      ! Inverse of the determinant
      DeterminantInverse = 1.D0 / Determinant

      ! Vertex position (solution of the system of equations)
      vCoordinates(1) = ( cSquaredRelativeDistance(iCand) * cRelativePosition(jCand,2) - cSquaredRelativeDistance(jCand) * &
      &                 cRelativePosition(iCand,2) ) * DeterminantInverse
      vCoordinates(2) = ( cSquaredRelativeDistance(jCand) * cRelativePosition(iCand,1) - cSquaredRelativeDistance(iCand) * &
      &                 cRelativePosition(jCand,1) ) * DeterminantInverse

      ! Initialization
      VertexValidity = .TRUE.

      ! Loop over the remaining candidates
      lCand = 1

      ! Plane inequality check
      DO WHILE( VertexValidity .AND. (lCand <= nCandidates) )

        ! Skip the current candidate
        IF( ( lCand /= iCand ) .AND. ( lCand /= jCand ) ) THEN
          ! Check if the vertex is closer to the central particle than to remaining candidates
          VertexValidity = ( DOT_PRODUCT( cRelativePosition(lCand,1:2), vCoordinates(1:2) ) <= cSquaredRelativeDistance(lCand) )
        END IF

        ! Next candidate
        lCand = lCand + 1

      END DO

      ! Add the vertex to the list of valid vertices
      IF( VertexValidity ) THEN

        ! Increment the number of valid vertices
        ValidVertices = ValidVertices + 1

        ! Check if the number of vertices exceeds the maximum allowed
        IF( ValidVertices > MaxVertices ) STOP "Program terminated: Maximum number of vertices exceeded!"

        ! Store the index of the candidate particles associated with the common vertex
        vIndex(ValidVertices,1) = iCand
        vIndex(ValidVertices,2) = jCand

        ! Convert the vertex position to the correct scale (this is a problem that arises from the use of relative positions)
        vRelativePosition(ValidVertices,1) = 0.5D0 * vCoordinates(1)
        vRelativePosition(ValidVertices,2) = 0.5D0 * vCoordinates(2)

      END IF

    END IF

  END DO

END DO

! Number of vertices
nVertices = ValidVertices

! Check if the number of vertices is less than 3 (minimum required for a polygon)
IF( nVertices < 3 ) THEN
  WRITE( *, "(3G0)" ) "Program terminated: Less than 3 vertices (", nVertices, ") found during the Voronoi construction!"
  STOP
END IF

! Identify neighboring points
cEdges = 0
DO iVertex = 1, nVertices
  ! Number of vertices per edge is the number of times a candidate was used to form an edge
  cEdges(vIndex(iVertex,1)) = cEdges(vIndex(iVertex,1)) + 1
  cEdges(vIndex(iVertex,2)) = cEdges(vIndex(iVertex,2)) + 1
END DO

! Calculate the number edges
VertexValidity = .TRUE.
nEdges = 0
DO iCand = 1, nCandidates
  ! Increment the number of edges
  IF( cEdges(iCand) > 0 ) THEN
    nEdges = nEdges + 1
    IF( cEdges(iCand) /= 2 ) THEN
      VertexValidity = .FALSE.
    END IF
  END IF
END DO

! Check inconsistencies
IF( .NOT. VertexValidity ) THEN
  WRITE( *, "(G0)" ) "Program terminated: Vertex degenaracy found during the Voronoi construction!"
  STOP
END IF

RETURN

END SUBROUTINE VoronoiConstruction2D

! ********************************************************************************************************************************* !
! Subroutine to calculate the Voronoi neighbors of a given particle in 3D.                                                          !
! ********************************************************************************************************************************* !
! See Allen and Tildesley, Computer Simulation of Liquids, 1987, for more information or visit                                      !
! <https://server.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/>.                                                      !
! ********************************************************************************************************************************* !
SUBROUTINE VoronoiConstruction3D( nCandidates, nVertices, nEdges, nFaces, cRelativePosition, cSquaredRelativeDistance, &
&                                 nFaceEdges, vRelativePosition, vIndex )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )                :: iCand, jCand, kCand ! Counters (candidates)
INTEGER( Kind= Int64 )                :: lCand               ! Counter (test candidates)
INTEGER( Kind= Int64 )                :: ValidVertices       ! Number of valid vertices
INTEGER( Kind= Int64 )                :: iVertex             ! Counters (vertices)
INTEGER( Kind= Int64 ), INTENT( IN )  :: nCandidates         ! Number of candidates
INTEGER( Kind= Int64 ), INTENT( OUT ) :: nVertices           ! Number of vertices
INTEGER( Kind= Int64 ), INTENT( OUT ) :: nEdges              ! Number of edges
INTEGER( Kind= Int64 ), INTENT( OUT ) :: nFaces              ! Number of faces

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), INTENT( OUT ), DIMENSION( : )    :: nFaceEdges ! Number of edges per face
INTEGER( Kind= Int64 ), INTENT( OUT ), DIMENSION( :, : ) :: vIndex     ! Index of the vertices

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Determinant        ! Determinant (check if the planes intersect)
REAL( Kind= Real64 ) :: DeterminantInverse ! Inverse of determinant

! REAL VARIABLES (SCALAR,PARAMETER)
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-10 ! Tolerance

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                   :: PlaneIntersect           ! Plane intersection coefficients
REAL( Kind= Real64 ), DIMENSION( 3 )                   :: NormalVector             ! Normal vector to the plane
REAL( Kind= Real64 ), DIMENSION( 3 )                   :: vCoordinates             ! Vertex coordinates
REAL( Kind= Real64 ), INTENT( IN ), DIMENSION( : )     :: cSquaredRelativeDistance ! Squared relative distance
REAL( Kind= Real64 ), INTENT( IN ), DIMENSION( :, : )  :: cRelativePosition        ! Relative position of the candidates
REAL( Kind= Real64 ), INTENT( OUT ), DIMENSION( :, : ) :: vRelativePosition        ! Relative position of the vertices

! LOGICAL VARIABLES
LOGICAL :: VertexValidity = .FALSE.

! Initialization
ValidVertices = 0

! Examination of each possible vertex defined by the intersection of 3 planes
DO iCand = 1, nCandidates - 2

  DO jCand = iCand + 1, nCandidates - 1

    ! Cross product of the relative position of the candidates to obtain the normal vector
    NormalVector(1) = cRelativePosition(iCand,2) * cRelativePosition(jCand,3) - cRelativePosition(iCand,3) * &
    &                 cRelativePosition(jCand,2)
    NormalVector(2) = cRelativePosition(iCand,3) * cRelativePosition(jCand,1) - cRelativePosition(iCand,1) * &
    &                 cRelativePosition(jCand,3)
    NormalVector(3) = cRelativePosition(iCand,1) * cRelativePosition(jCand,2) - cRelativePosition(iCand,2) * &
    &                 cRelativePosition(jCand,1)

    ! Plane intersection coefficients
    PlaneIntersect(:) = DOT_PRODUCT( cRelativePosition(jCand,1:3), cRelativePosition(jCand,1:3) ) * cRelativePosition(iCand,:) - &
    &                   DOT_PRODUCT( cRelativePosition(iCand,1:3), cRelativePosition(iCand,1:3) ) * cRelativePosition(jCand,:) 

    DO kCand = jCand + 1, nCandidates

      ! Calculate the determinant and check if the planes intersect
      Determinant = DOT_PRODUCT( NormalVector(1:3), cRelativePosition(kCand,1:3) )

      ! Condition for plane intersection (non-zero determinant means unique solution or a unique point instead of a line)
      IF( DABS( Determinant ) > Tolerance ) THEN

        ! Inverse of the determinant
        DeterminantInverse = 1.D0 / Determinant

        ! Vertex position (solution of the system of equations)
        vCoordinates(1) = ( cSquaredRelativeDistance(kCand) * NormalVector(1) + cRelativePosition(kCand,2) * PlaneIntersect(3) - &
        &                 cRelativePosition(kCand,3) * PlaneIntersect(2) ) * DeterminantInverse
        vCoordinates(2) = ( cSquaredRelativeDistance(kCand) * NormalVector(2) + cRelativePosition(kCand,3) * PlaneIntersect(1) - &
        &                 cRelativePosition(kCand,1) * PlaneIntersect(3) ) * DeterminantInverse
        vCoordinates(3) = ( cSquaredRelativeDistance(kCand) * NormalVector(3) + cRelativePosition(kCand,1) * PlaneIntersect(2) - &
        &                 cRelativePosition(kCand,2) * PlaneIntersect(1) ) * DeterminantInverse

        ! Initialization
        VertexValidity = .TRUE.

        ! Loop over the remaining candidates
        lCand = 1

        ! Plane inequality check
        DO WHILE( VertexValidity .AND. (lCand <= nCandidates) )

          ! Skip the current candidate
          IF( ( lCand /= iCand ) .AND. ( lCand /= jCand ) .AND. ( lCand /= kCand ) ) THEN
            ! Check if the vertex is closer to the central particle than to remaining candidates
            VertexValidity = ( DOT_PRODUCT( cRelativePosition(lCand,1:3), vCoordinates(1:3) ) <= cSquaredRelativeDistance(lCand) )
          END IF

          ! Next candidate
          lCand = lCand + 1

        END DO

        ! Add the vertex to the list of valid vertices
        IF( VertexValidity ) THEN

          ! Increment the number of valid vertices
          ValidVertices = ValidVertices + 1

          ! Check if the number of vertices exceeds the maximum allowed
          IF( ValidVertices > MaxVertices ) STOP "Program terminated: Maximum number of vertices exceeded!"

          ! Store the index of the candidate particles associated with the common vertex
          vIndex(ValidVertices,1) = iCand
          vIndex(ValidVertices,2) = jCand
          vIndex(ValidVertices,3) = kCand

          ! Convert the vertex position to the correct scale (this is a problem that arises from the use of relative positions)
          vRelativePosition(ValidVertices,1) = 0.5D0 * vCoordinates(1)
          vRelativePosition(ValidVertices,2) = 0.5D0 * vCoordinates(2)
          vRelativePosition(ValidVertices,3) = 0.5D0 * vCoordinates(3)

        END IF

      END IF

    END DO

  END DO

END DO

! Number of vertices
nVertices = ValidVertices

! Check if the number of vertices is less than 4 (minimum required for a tetrahedron)
IF( nVertices < 4 ) THEN
  WRITE( *, "(3G0)" ) "Program terminated: Less than 4 vertices (", nVertices, ") found during the Voronoi construction!"
  STOP
END IF

! Identify neighboring points
nFaceEdges = 0
DO iVertex = 1, nVertices
  ! Number of edges per face is the number of times a candidate was used to form a face (e.g. if the candidate was considered 6 times, than the face it shares with the central particle has 6 vertices or 6 edges)
  nFaceEdges(vIndex(iVertex,1)) = nFaceEdges(vIndex(iVertex,1)) + 1
  nFaceEdges(vIndex(iVertex,2)) = nFaceEdges(vIndex(iVertex,2)) + 1
  nFaceEdges(vIndex(iVertex,3)) = nFaceEdges(vIndex(iVertex,3)) + 1
END DO

! Calculate the number of faces and edges
nFaces = 0
nEdges = 0
DO iCand = 1, nCandidates
  ! Increment the number of faces
  IF( nFaceEdges(iCand) > 0 ) nFaces = nFaces + 1
  ! Increment the number of edges
  nEdges = nEdges + nFaceEdges(iCand)
END DO

! Check inconsistencies (an edge is shared by two faces)
IF( MOD ( nEdges, 2 ) /= 0 ) THEN
  WRITE( *, "(3G0)" ) "Program terminated: Non-integer number of edges (", nEdges, ") found during the Voronoi construction!"
  STOP
END IF

! Eliminate the double counting of edges
nEdges = nEdges / 2

! Check the Euler's relation (V - E + F = 2)
IF( ( nVertices - nEdges + nFaces ) /= 2 ) THEN
  WRITE( *, "(7G0)" ) "Program terminated: Euler's relation (", nVertices, " - ", nEdges, " + ", nFaces, " = 2) not satisfied!"
  STOP
END IF

RETURN

END SUBROUTINE VoronoiConstruction3D

END MODULE Neighbors