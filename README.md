# Hex1Mesh
Hexahedral Mesh repository. Includes generation, visualization, basefunctions, CutFEM suites, etc.

## Functions
**Hex1Mesh** includes the following functions:

- Hex1Mesh() - The constructor, creates the Hex1Mesh object.
- Neighbors() - Returns a m-by-4 matrix containing element indices to the neighbors.
- vizMesh()
- ElementVolume()
- RefineLocal()
- CutP1()
- vizP1Surf()
- TriangulateP1()
- CutP2()
- baseHexP1()

----------

**TODO:**

- baseHexQuad()


## Mesh numbering
The mesh is numbered according to the figure below. The red and cyan colors are used to explain the refinement algorithm.
![](http://www.mirzacenanovic.com/wp-content/uploads/2015/01/2015-01-09-14_50_26-Figure-1_-Hex1Mesh.png)

