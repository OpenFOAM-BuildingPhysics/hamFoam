Info<< "Reading buildingMaterial Information\n" << endl;

labelList faceType(mesh.nFaces(), -1);
labelList cellType(mesh.nCells(), -1);
//const labelList& own = mesh.faceOwner();

PtrList<dictionary> Materials(transportProperties.lookup("buildingMaterials"));
