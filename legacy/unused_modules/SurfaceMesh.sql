
DROP TABLE #outer_surface_vertices
DROP TABLE #outer_surface_edges
DROP TABLE #outer_surface_triangles

SELECT VertexID INTO #outer_surface_vertices
FROM TFaceVertices
WHERE VertexID IN (
  SELECT VertexID
  FROM TFaceVertices
  WHERE FaceID IN (
    SELECT FaceID FROM OuterTFaces
  )
)

SELECT DISTINCT EdgeID, v0, v1 INTO #outer_surface_edges
FROM Edges
  JOIN TFaceVertices AS A ON v0 = A.VertexID
  JOIN TFaceVertices AS B ON v1 = B.VertexID
  WHERE A.FaceID = B.FaceID
    AND A.FaceID IN (
      SELECT FaceID FROM OuterTFaces
    )

SELECT VertexID, x, y, z FROM Vertices WHERE VertexID IN (SELECT VertexID FROM #outer_surface_vertices)
UNION ALL
SELECT NodeID, x, y, z FROM Nodes
  WHERE NodeID IN (
    SELECT EdgeID FROM #outer_surface_edges
  )

SELECT TriID INTO #outer_surface_triangles
FROM TFaceTriangles
WHERE FaceID IN (
    SELECT FaceID FROM OuterTFaces
)

SELECT TriID, D.v0, D.v1, D.v2, A.EdgeID AS v3, B.EdgeID AS v4, C.EdgeID AS v5
FROM (
SELECT TriID, v0, v1, v2,
  v3 = CASE WHEN v0 < v1 THEN v0 ELSE v1 END,
  v4 = CASE WHEN v0 < v1 THEN v1 ELSE v0 END,
  v5 = CASE WHEN v1 < v2 THEN v1 ELSE v2 END,
  v6 = CASE WHEN v1 < v2 THEN v2 ELSE v1 END,
  v7 = CASE WHEN v0 < v2 THEN v0 ELSE v2 END,
  v8 = CASE WHEN v0 < v2 THEN v2 ELSE v0 END
FROM TFaceTriangles
WHERE TriID IN (SELECT TriID FROM #outer_surface_triangles)
) AS D
  JOIN Edges AS A ON A.v0 = D.v3 AND A.v1 = D.v4
  JOIN Edges AS B ON B.v0 = D.v5 AND B.v1 = D.v6
  JOIN Edges AS C ON C.v0 = D.v7 AND C.v1 = D.v8
  WHERE A.EdgeID IN (SELECT EdgeID FROM #outer_surface_edges)
    AND B.EdgeID IN (SELECT EdgeID FROM #outer_surface_edges)
    AND C.EdgeID IN (SELECT EdgeID FROM #outer_surface_edges)

select max(x) as "max_x", max(y) as "max_y", max(z) as "max_z", min(x) as "min_x", min(y) as "min_y", min(z) as "min_z"
from (SELECT VertexID, x, y, z FROM Vertices WHERE VertexID IN (SELECT VertexID FROM #outer_surface_vertices)
UNION ALL
SELECT NodeID, x, y, z FROM Nodes
  WHERE NodeID IN (
    SELECT EdgeID FROM #outer_surface_edges
	)
)

select max(x), max(y), max(z), min(x), min(y), min(z)
from TVertices