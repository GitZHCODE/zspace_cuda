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

#include <headers/zCudaToolsets/energy/zTsSolarOcclusionGPU.h>
using namespace std;
#include <chrono>

namespace zSpace
{

	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSolarOcclusionGPU::zTsSolarOcclusionGPU() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSolarOcclusionGPU::~zTsSolarOcclusionGPU() {}

	//---- SET METHODS

	ZSPACE_INLINE void zTsSolarOcclusionGPU::setBuffer(float* _bPositions, float* _bVals, int _numBuffer)
	{
		bufferPositions = new float[_numBuffer];;
		std::copy(_bPositions, _bPositions + (_numBuffer), bufferPositions);

		bufferValues = new float[_numBuffer];
		std::copy(_bVals, _bVals + (_numBuffer), bufferValues);

		num_Buffer = _numBuffer;
	}
		
	ZSPACE_INLINE void zTsSolarOcclusionGPU::setFacePositions(float* _fPositions, int _numFace)
	{
		facePositions = new float[_numFace];;
		std::copy(_fPositions, _fPositions + _numFace, facePositions);

		numFaces = _numFace;

		setMemory();

		std::copy(bufferPositions, bufferPositions + (num_Buffer), buffer_facePositions);
		std::copy(facePositions, facePositions + numFaces, buffer_facePositions + num_Buffer);	
	}

	//---- GET METHODS
	ZSPACE_INLINE int zTsSolarOcclusionGPU::numBuffer()
	{
		return num_Buffer;
	}

	ZSPACE_INLINE int zTsSolarOcclusionGPU::numFacePositions()
	{
		return numFaces;
	}

	ZSPACE_INLINE float* zTsSolarOcclusionGPU::getRawBuffer()
	{
		return bufferPositions;
	}

	ZSPACE_INLINE float* zTsSolarOcclusionGPU::getRawBufferValues()
	{
		return bufferValues;
	}

	ZSPACE_INLINE float* zTsSolarOcclusionGPU::getRawFacePositions()
	{
		return facePositions;
	}

	ZSPACE_INLINE float* zTsSolarOcclusionGPU::getRawBuffer_FacePositions()
	{
		return buffer_facePositions;
	}



	//---- COMPUTE METHODS

	ZSPACE_INLINE void zTsSolarOcclusionGPU::setMemory()
	{
		if ((num_Buffer + numFaces) < memSize) return;
		else
		{
			while (memSize < (num_Buffer + numFaces)) memSize += d_MEMORYMULTIPLIER;
		
			buffer_facePositions = new float[memSize];
		
			// set to  Num Normals
			//colors = new float[memSize];
		}
	}
}
