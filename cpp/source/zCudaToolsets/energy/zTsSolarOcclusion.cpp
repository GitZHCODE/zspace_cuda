#include <headers/zCudaToolsets/energy/zTsSolarOcclusion.h>

namespace zSpace
{

	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSolarOcclusion::zTsSolarOcclusion() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSolarOcclusion::~zTsSolarOcclusion() {}

	//---- SET METHODS
	ZSPACE_INLINE void zTsSolarOcclusion::setSunVec(zVector &_sunVec)
	{
		sunVector = &_sunVec;
	}

	ZSPACE_INLINE void zTsSolarOcclusion::setSunCenter(zVector &_sunCenter)
	{
		sunCenter = &_sunCenter;
	}

	ZSPACE_INLINE void zTsSolarOcclusion::setResolution(int & _resolutionX, int& _resolutionY)
	{
		resolutionX = _resolutionX;
		resolutionY = _resolutionY;
	}

	ZSPACE_INLINE void zTsSolarOcclusion::setBuildingMesh(zObjMesh _buildingMesh)
	{
		buildingMesh = _buildingMesh;
	}


	//---- GET METHODS
	ZSPACE_INLINE zBoolArray zTsSolarOcclusion::getOcclusion()
	{
		return vertexOcclusion;
	}

	ZSPACE_INLINE zIntArray zTsSolarOcclusion::getOcclusionValue()
	{
		return occlusionVal;
	}

	ZSPACE_INLINE zObjMeshScalarField zTsSolarOcclusion::getField()
	{
		return oMesh_Field;
	}

	ZSPACE_INLINE zTransform zTsSolarOcclusion::getTransform()
	{
		return transform;
	}

	ZSPACE_INLINE zTransform zTsSolarOcclusion::getTransformLocal()
	{
		return transformLocal;
	}



	//---- COMPUTE METHODS
	 
	//global
	ZSPACE_INLINE void zTsSolarOcclusion::init()
	{
		zFnMesh fnMesh(buildingMesh);
		numVertices = fnMesh.numVertices();
		numFaces = fnMesh.numPolygons();

		occlusionVal.assign(numVertices, 0);
	}

	ZSPACE_INLINE void zTsSolarOcclusion::computeOcclusion()
	{
		init();
		sortFaceByDist();
		computeTransform();
		flatternBuilding();
		computeFieldBBox();
		createField();

		for (int i = 0; i < numFaces; i++)
		{
			computeVertexValue(i);
			updateField(i);
		}

		zFnMeshScalarField f(oMesh_Field);
		f.updateColors();
	}


	//sun
	ZSPACE_INLINE void zTsSolarOcclusion::computeTransform()
	{
		zVector O = *sunCenter;
		zVector Z = *sunVector;
		zTransform t;

		Z.normalize();
		zVector basis(0, 1, 0);
		zVector X = basis ^ Z;

		zVector Y = Z ^ X;
		Y.normalize();

		X = Y ^ Z;
		X.normalize();

		t.setIdentity();
		t(0, 0) = X.x; t(0, 1) = X.y; t(0, 2) = X.z;
		t(1, 0) = Y.x; t(1, 1) = Y.y; t(1, 2) = Y.z;
		t(2, 0) = Z.x; t(2, 1) = Z.y; t(2, 2) = Z.z;
		t(3, 0) = O.x; t(3, 1) = O.y; t(3, 2) = O.z;

		//output
		transform = t;
		transformLocal.setIdentity();
	}

	//field
	ZSPACE_INLINE void zTsSolarOcclusion::computeFieldBBox()
	{
		zFnMesh fnMesh(buildingMesh);
		fnMesh.getBounds(fieldMinBB, fieldMaxBB);
	}

	ZSPACE_INLINE void zTsSolarOcclusion::createField()
	{
		//get field
		zFnMeshScalarField fnField(oMesh_Field);

		fieldMinBB += zVector(-1, -1, 0);
		fieldMaxBB += zVector(1, 1, 0);

		int fieldSize = resolutionX * resolutionY;
		unitX = fieldMaxBB.x / resolutionX;
		unitY = fieldMaxBB.y / resolutionY;
		fnField.create(fieldMinBB, fieldMaxBB, resolutionX, resolutionY, true, false);

		//init field value
		zFloatArray value;
		value.assign(fieldSize, 0);
		fnField.setFieldValues(value);
	}

	ZSPACE_INLINE void zTsSolarOcclusion::updateField(int currentfaceID)
	{
		//get min bounds from projected face 2d
		zPoint minBB, maxBB;
		zFnMesh fnMesh(buildingMesh);
		zItMeshFace face(buildingMesh, faceID[currentfaceID]);
		zPointArray poly;

		face.getVertexPositions(poly);
		core.getBounds(poly, minBB, maxBB);

		zFnMeshScalarField sf(oMesh_Field);

		double fieldUnitX, fieldUnitY;
		int X, Y;
		sf.getUnitDistances(fieldUnitX, fieldUnitY);
		sf.getResolution(X, Y);

		int _index_X = floor((minBB.x - fieldMinBB.x) / fieldUnitX);
		int _index_Y = floor((minBB.y - fieldMinBB.y) / fieldUnitY);

		int _index = _index_Y * X + _index_X;

		zItMeshScalarField s_min(oMesh_Field, minBB, false);
		int minID = s_min.getId();
		int minID_x, minID_y;
		s_min.getIndices(minID_x, minID_y);
		
		zItMeshScalarField s_max(oMesh_Field, maxBB, false);
		int maxID = s_max.getId();
		int maxID_x, maxID_y;
		s_max.getIndices(maxID_x, maxID_y);

		//check in poly and set value
		for (int i = minID_x; i <= maxID_x; i++)
		{
			for (int j = minID_y; j <= maxID_y; j++)
			{
				zItMeshScalarField s(oMesh_Field, i, j);
				zPoint p = s.getPosition();
				bool IsInPoly = core.pointInPlanarPolygon(p, poly,zVector(0,0,1));
				if (IsInPoly)
				{
					s.setValue(1);
				}
			}
		}
	}

	//building
	ZSPACE_INLINE void zTsSolarOcclusion::sortFaceByDist()
	{
		zDoubleArray dist(numFaces);
		zVector O = *sunCenter;
		zVector Z = *sunVector;
		zFnMesh fnMesh(buildingMesh);
		zPointArray centers;
		fnMesh.getCenters(zFaceData, centers);

		//face center distance to sunPlane
		for (size_t i = 0; i < numFaces; i++)
		{
			dist[i] = core.minDist_Point_Plane(centers[i], O, Z);
		}

		//std sort faceId by distance
		vector<size_t> idx(numFaces);
		iota(idx.begin(), idx.end(), 0);
		stable_sort(idx.begin(), idx.end(),
			[&dist](size_t a, size_t b) {return dist[a] < dist[b]; });

		//output
		faceID = idx;
	}

	ZSPACE_INLINE void zTsSolarOcclusion::flatternBuilding()
	{
		//project buildingMesh
		zVector O = *sunCenter;
		zVector Z = *sunVector;
		zFnMesh fnMesh(buildingMesh);

		for (int i = 0; i < numVertices; i++)
		{
			double dist = core.minDist_Point_Plane(fnMesh.getRawVertexPositions()[i], O, Z);
			zVector target = Z * dist * -1;
			fnMesh.getRawVertexPositions()[i] += target;
		}

		//transform buildingMesh
		fnMesh.setTransform(transform, true, false);
		fnMesh.setTransform(transformLocal, true, true);
	}

	ZSPACE_INLINE void zTsSolarOcclusion::computeVertexValue(int currentfaceID)
	{
		zItMeshFace face(buildingMesh, faceID[currentfaceID]);
		zPointArray v; 
		face.getVertexPositions(v);

		zIntArray pId;
		face.getVertices(pId);

		for (int i = 0; i < v.size(); i++)
		{
			zItMeshScalarField s(oMesh_Field, v[i], false);
			int value = s.getValue();

			if (!value)
			{
				occlusionVal[pId[i]]++;
			}
		}
	}

}
