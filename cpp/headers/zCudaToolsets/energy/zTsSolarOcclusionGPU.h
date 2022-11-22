// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author :  Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Taizhong Chen <Taizhong.Chen@zaha-hadid.com>
//

#ifndef ZSPACE_TS_SOLAR_OCCLUSION_GPU
#define ZSPACE_TS_SOLAR_OCCLUSION_GPU

#pragma once
#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zDomain.h>
#include<headers/zCore/utilities/zUtilsCore.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

namespace zSpace
{
	//struct zNorm_SunVec
	//{
	//	zVector norm;
	//	zVector sunVec;
	//};

	class ZSPACE_CUDA_TOOLSET zTsSolarOcclusionGPU

	{
	protected:
		/*!	\brief pointer to container of buffer positions	*/
		float* bufferPositions;
		
		/*!	\brief pointer to container of buffer values	*/
		int* bufferValues;

		/*!	\brief number of buffer positions in the container, num_Buffer = resX by resY by 3 */
		int num_Buffer;

		/*!	\brief pointer to container of building quad face connected vertices and face distance	*/
		float* facePositions;

		/*!	\brief number of building faces in the container, numFaces = numFace by 4 by 4	*/
		int numFaces;

		/*!	\brief pointer to container of occlusion values */
		float* occlusionVal;

		/*!	\brief pointer to container of normals + building quad face connected vertices and face distance */
		float* buffer_facePositions;

		/*!	\brief size of container for normals, building quad face connected vertices, face distance */
		int memSize;

		//	//--------------------------
		//	//---- COMPUTE METHODS
		//	//--------------------------

	public:

		zUtilsCore core;

		//	/*!	\brief core utilities Object  */
			//--------------------------
			//---- CONSTRUCTOR
			//--------------------------

			/*! \brief Default constructor.
			*	\since version 0.0.4
			*/
		zTsSolarOcclusionGPU();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		~zTsSolarOcclusionGPU();

		//	//--------------------------
		//	//---- SET METHODS
		//	//--------------------------	
		//		   
		ZSPACE_CUDA_CALLABLE_HOST void setBuffer(float* _bPositions, int* _bVals, int _numBuffer);

		ZSPACE_CUDA_CALLABLE_HOST void setFacePositions(float* _fPositions,  int _numFace);

		//	//--------------------------
		//	//---- GET METHODS
		//	//--------------------------	

		ZSPACE_CUDA_CALLABLE int numBuffer();

		ZSPACE_CUDA_CALLABLE int numFacePositions();

		ZSPACE_CUDA_CALLABLE float* getRawBuffer();

		ZSPACE_CUDA_CALLABLE int* getRawBufferValues();

		ZSPACE_CUDA_CALLABLE float* getRawFacePositions();

		ZSPACE_CUDA_CALLABLE float* getRawBuffer_FacePositions();
		

		//	//--------------------------
		//	//---- COMPUTE METHODS
		//	//--------------------------
		

	protected:

		ZSPACE_CUDA_CALLABLE_HOST void setMemory();

	
	};
}

#if defined(ZSPACE_CUDA_STATIC_LIBRARY)  || defined(ZSPACE_CUDA_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/energy/zTsSolarOcclusionGPU.cpp>
#endif

#endif