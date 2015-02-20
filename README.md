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
- vizP2Surf() (includes P2 triangulation)
- boundaryInds()
- baseHexP1()

## Demo
    % Create mesh
	x0 = -1;
	x1 = 1;
	y0 = -1;
	y1 = 1;
	z0 = -1;
	z1 = 1;
	ne = 2;
	nxe = ne;
	nye = ne;
	nze = ne;
	H = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
	H.vizMesh('ElementNumbers','NodeNumbers');
	
	% Create surface
	R = 0.89;
	xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
	surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;
	phi = surfaceFunction(H.XC, H.YC, H.ZC);
	
	H.CutP2(phi,0)
	hv = H.vizP2Surf('FaceColor','c','EdgeColor','k','nP2Ele',5); axis tight
	

## Mesh numbering
The mesh is numbered according to the figure below. The red and cyan colors are used to explain the refinement algorithm.
![](http://www.mirzacenanovic.com/wp-content/uploads/2015/01/2015-01-09-14_50_26-Figure-1_-Hex1Mesh.png)

### Examples
- P2 surface:
![](https://raw.githubusercontent.com/cenmir/Hex1Mesh/master/Hex1Mesh/Examples/P2.png)

- Issue with mid nodes not being in-plane:
![](https://raw.githubusercontent.com/cenmir/Hex1Mesh/master/Hex1Mesh/Examples/noInPlaneSearch.png)

- In plane projection:
![](https://raw.githubusercontent.com/cenmir/Hex1Mesh/master/Hex1Mesh/Examples/Tangential%20Projection.png)
(Blue is original search direction, black is in-plane search direction, green is face normal)