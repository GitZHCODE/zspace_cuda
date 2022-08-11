#include <headers/zCudaToolsets/field/zTsPixelField.h>




//-----
// zTsMeshField
//-----
namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsPixelField::zTsPixelField() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsPixelField::~zTsPixelField() {}

	//---- SET METHODS


	ZSPACE_INLINE void zTsPixelField::setPixels(const float* _vPositions, const float* _vColors, float* _pixelValues, const int* _polyconnects, int _numVertices, int _numPolygons)
	{
		numV = _numVertices;
		vPositions = new float[_numVertices * 3];;
		std::copy(_vPositions, _vPositions + (_numVertices * 3), vPositions);

		pixelValues = new float[_numVertices];;
		std::copy(_pixelValues, _pixelValues + (_numVertices * 1), pixelValues);

		numP = _numPolygons;
		polyConnects = new int[_numPolygons * 4];
		std::copy(_polyconnects, _polyconnects + (_numPolygons * 4), polyConnects);

		numVoxs = _numVoxels;
		setMemory(numVoxs * 4 + MAX_POINTCLOUD);


		int numPositions = numVoxs * 3;
		std::copy(voxelPositions, voxelPositions + numPositions, voxels);
		std::copy(voxelData, voxelData + numVoxs, voxels + (numVoxs * 3));
	}

	ZSPACE_INLINE void zTsVoxelField::setPointcloud(float* _pointCloud, int _numPoints)
	{
		pointCloud = new float[_numPoints];
		std::copy(_pointCloud, _pointCloud + (_numPoints * 1), pointCloud);
		numPoints = _numPoints;

		if (numPoints > MAX_POINTCLOUD)
		{
			setMemory((numVoxs * 4) + numPoints);

			std::copy(voxelPositions, voxelPositions + (numVoxs * 3), voxels);
			std::copy(voxelData, voxelData + numVoxs, voxels + (numVoxs * 3));
		}

		
		std::copy(pointCloud, pointCloud + numPoints, voxels + (numVoxs * 4));
	}

	ZSPACE_INLINE void zTsVoxelField::setDomain_Colors(zDomainColor& _dColor)
	{
		dColor = _dColor;
	}


	//---- GET METHODS

	ZSPACE_INLINE int zTsVoxelField::numVoxels()
	{
		return numVoxs;
	}

	ZSPACE_INLINE int zTsVoxelField::numCloudPoints()
	{
		return numPoints;
	}

	ZSPACE_INLINE int zTsVoxelField::getMemSize()
	{
		return memSize;
	}

	ZSPACE_INLINE float* zTsVoxelField::getRawVoxels()
	{
		return voxels;
	}

	ZSPACE_INLINE float* zTsVoxelField::getRawColors()
	{
		return voxelColors;
	}

	ZSPACE_INLINE zDomainColor zTsVoxelField::getDomain_Colors()
	{
		return dColor;
	}


	
	//---- COMPUTE METHODS

	ZSPACE_INLINE void zTsVoxelField::computeBin()
	{
		
	}


	//---- PROTECTED METHODS


	ZSPACE_INLINE void zTsVoxelField::setMemory(int nSize)
	{
		if ((nSize) < memSize) return;
		else
		{
			while (memSize < (nSize)) memSize += d_MEMORYMULTIPLIER;
			voxels = new float[memSize];	

			voxelColors = new float[memSize];
			
			printf("\n cpu  MemSize %i \n", memSize);
		}
	}


}
