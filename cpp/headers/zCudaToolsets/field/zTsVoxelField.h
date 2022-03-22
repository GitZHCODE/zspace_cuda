// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Federico Borello <federico.borello@zaha-hadid.com>
//

#ifndef ZSPACE_TS_VOXEL_FIELD
#define ZSPACE_TS_VOXEL_FIELD

#pragma once

#define MAX_POINTCLOUD 512 * 50

#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zDomain.h>
//#include<headers/zCore/utilities/zUtilsCore.h>

namespace zSpace
{
	

	class ZSPACE_CUDA zCdGraph
	{
		

	public:

		/*!	\brief pointer to container of vertex positions	*/
		float* vertexPositions;

		/*!	\brief pointer to container of vertex positions	*/
		int* edgeConnects;

		/*!	\brief number of vertices in the container	*/
		int nV;

		/*!	\brief number of vertices in the container	*/
		int nE;


		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST zCdGraph();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST ~zCdGraph();

		//--------------------------
		//---- SET METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE_HOST void setGraph(const float* _vertexPositions, const int* _edgeConnects, int _numVertices, int _numEdges);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		ZSPACE_CUDA_CALLABLE int numVertices();

		ZSPACE_CUDA_CALLABLE int numEdges();

		ZSPACE_CUDA_CALLABLE float* getRawVertexPositions();

		ZSPACE_CUDA_CALLABLE int* getRawEdgeConnects();

	};

	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsSolarAnalysis
	*	\brief A tool set to do solar analysis.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_CUDA zTsVoxelField
	{
	private:



		/*!	\brief pointer to container of voxel positions + data + pointcloud	*/
		float* voxels;

		/*!	\brief pointer to container of voxel positions	*/
		float* voxelPositions;

		/*!	\brief pointer to container of voxel colors*/
		float* voxelColors;

		/*!	\brief pointer to container of pont cloud*/
		float* pointCloud;

		/*!	\brief pointer to container of voxel data	*/
		float* voxelData;

		/*!	\brief number of voxels in the container	*/
		int numVoxs;

		/*!	\brief number of points in the container	*/
		int numPoints;

		/*!	\brief color domain	*/
		zDomainColor dColor;



		/*!	\brief memory size of container */
		int memSize;



	public:

		/*!	\brief core utilities Object  */
		//zUtilsCore coreUtils;


		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST zTsVoxelField();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST ~zTsVoxelField();

		//--------------------------
		//---- SET METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE_HOST void setVoxels(float* _voxelPositions, float* _voxelData, int _numVoxels);

		ZSPACE_CUDA_CALLABLE_HOST void setPointcloud(float* _pointCloud, int _numPoints);

		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Colors(zDomainColor& _dColor);

		//--------------------------
		//---- GET METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE int numVoxels();

		ZSPACE_CUDA_CALLABLE int getMemSize();

		ZSPACE_CUDA_CALLABLE int numCloudPoints();

		ZSPACE_CUDA_CALLABLE float* getRawVoxels();

		ZSPACE_CUDA_CALLABLE float* getRawColors();

		ZSPACE_CUDA_CALLABLE zDomainColor getDomain_Colors();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE void computeBin();
				

		//--------------------------
		//---- DISPLAY METHODS
		//--------------------------	

		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------	
	protected:

		ZSPACE_CUDA_CALLABLE_HOST void setMemory( int nSize);

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/field/zTsVoxelField.cpp>
#endif

#endif