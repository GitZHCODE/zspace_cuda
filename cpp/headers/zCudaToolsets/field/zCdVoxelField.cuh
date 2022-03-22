// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//

#ifndef ZSPACE_CD_VOXELFIELD
#define ZSPACE_CD_VOXELFIELD

#pragma once




#include<headers/zCudaToolsets/base/zCdUtilities.cuh>
#include<headers/zCudaToolsets/field/zTsVoxelField.h>

using namespace zSpace;

//---- DEVICE VARIABLES

float* d_voxels;
float* d_colors;

int d_MemSize;

//---- CUDA HOST DEVICE METHODS 

ZSPACE_CUDA_CALLABLE float ofMap(float value, float inputMin, float inputMax, float outputMin, float outputMax)
{
	return ((value - inputMin) / (inputMax - inputMin) * (outputMax - outputMin) + outputMin);
}

ZSPACE_CUDA_CALLABLE float minDist_Edge_Point(zVector& pt, zVector& e0, zVector& e1, zVector& closest_Pt)
{
	zVector n = (e1 - e0) ^ (zVector(0, 0, 1));
	n.normalize();
	closest_Pt = n * ((e0 - pt) * n);
	closest_Pt += pt;

	float len = e0.distanceTo(e1);

	zVector ed = (e0 - e1) / len;
	double param = (closest_Pt - e1) * ed;



	if (param > len)param = len;
	if (param < 0) param = 0;

	closest_Pt = e1 + ed * param;

	return closest_Pt.distanceTo(pt);
}


ZSPACE_EXTERN void cleanDeviceMemory()
{
	// Free memory.
	cudaFree(d_voxels);
	cudaFree(d_colors);

}

ZSPACE_CUDA_CALLABLE_HOST void setDeviceMemory(int _nSize)
{
	if (_nSize < d_MemSize) return;
	else
	{
		while (d_MemSize < _nSize) d_MemSize += d_MEMORYMULTIPLIER;

		cleanDeviceMemory();

		checkCudaErrors(cudaMalloc((void**)&d_voxels, d_MemSize * FloatSize));
				
		checkCudaErrors(cudaMalloc((void**)&d_colors, d_MemSize * FloatSize));
	}
}

//---- CUDA KERNEL 

ZSPACE_CUDA_GLOBAL void computeFieldGraph_kernel(float* outField, float* outColors, int numVoxels, int numPoints)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < (numVoxels * 3) && i % 3 == 0)
	{
		zVector voxPos(outField[i + 0], outField[i + 1], outField[i + 2]);

		float d = 0.0;
		float tempDist = 100000;

		uint pointcloud_offset = (numVoxels * 4) - i;

		
		for (int o = i; o < i + numPoints; o += 6)
		{
			int j = o + pointcloud_offset;

			zVector p0(outField[j + 0], outField[j + 1], outField[j + 2]);
			zVector p1(outField[j + 3], outField[j + 4], outField[j + 5]);

			float edgeLen = p0.distanceTo(p1);
			if (edgeLen < EPS) continue;

			zVector closestPt;
			float r = minDist_Edge_Point(voxPos, p0, p1, closestPt);

			if (r < tempDist)
			{
				d = r;
				tempDist = r;
			}
		}		

		

		uint fieldvalue_offset = (numVoxels * 3) + i; /* (i == 0) ? (numVoxels) - i : (numVoxels) - i + floorf(i / 3)*/;
		outField[fieldvalue_offset] = d;

		//if (i == 0 || i == 3) printf("\n %1.2f | %1.2f %1.2f %1.2f", d, voxPos.x, voxPos.y, voxPos.z);

		// SDFs
		if (d < -1)
		{
			outColors[i + 0] = 0.25;
			outColors[i + 1] = 0.25;
			outColors[i + 2] = 0.25;

			//alpha
			outColors[fieldvalue_offset] = 0.25;
		}
		else if (d > 1)
		{

			outColors[i + 0] = 0.25;
			outColors[i + 1] = 0.25;
			outColors[i + 2] = 0.25;	

			//alpha
			outColors[fieldvalue_offset] = 0.25;
		}
		else
		{
			outColors[i + 0] = 0.25;
			outColors[i + 1] = 0.25;
			outColors[i + 2] = 0.25;

			//alpha
			outColors[fieldvalue_offset] = 0.25;
		}


	}

	
}

ZSPACE_CUDA_GLOBAL void computeFieldBlend(float* inField1, float* inField2, float* outField, int numFields, int numVoxels, float weight)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (i < numVoxels && i % 3 == 0)
	{
		uint fieldvalue_offset = (i == 0) ? numVoxels - i : numVoxels - i + floorf(i / 3);
		
		outField[fieldvalue_offset] = ((1 - weight) * inField1[fieldvalue_offset]) + (weight * inField2[fieldvalue_offset]);
	}

}

//---- launch KERNEL METHODS

ZSPACE_EXTERN bool cdpGraphBlends(zTsVoxelField& voxField)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_VOXELS = voxField.numVoxels();
	int NUM_POINTS = voxField.numCloudPoints();
	int h_MEMSIZE = voxField.getMemSize();
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//zDomainColor dColor = voxField.getDomain_Colors();
	setDeviceMemory(h_MEMSIZE);

	cudaEventRecord(start);

	// transfer memory to device

	checkCudaErrors(cudaMemcpy(d_voxels, voxField.getRawVoxels(), d_MemSize * FloatSize, cudaMemcpyHostToDevice));

	//checkCudaErrors(cudaMemcpy(d_colors, voxField.getRawColors(), d_MemSize * FloatSize, cudaMemcpyHostToDevice));

	// Launch Kernel
	
	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_MemSize / (double)block.x));
	
	computeFieldGraph_kernel << < grid, block >> > (d_voxels, d_colors, NUM_VOXELS, NUM_POINTS);

	checkCudaErrors(cudaGetLastError());


	// transfer memory to host
	checkCudaErrors(cudaMemcpy(voxField.getRawColors(), d_colors, d_MemSize * FloatSize, cudaMemcpyDeviceToHost));

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("\n gpu %1.8f ms \n", milliseconds);

	printf("\n gpu  d_MemSize %i \n", d_MemSize);

	return true;


}


#endif

// ZERO COPY METHOD
// build map example 
//cudaSetDeviceFlags(cudaDeviceMapHost);
//cudaHostAlloc(&sAnalysis.norm_sunVecs, NUM_ANGLES * sizeof(zNorm_SunVec), cudaHostAllocMapped);
//cudaHostGetDevicePointer(&d_norm_sunVecs, sAnalysis.norm_sunVecs, 0);

//cudaSetDeviceFlags(cudaDeviceMapHost);
//cudaHostAlloc(&sAnalysis.solarAngles, NUM_ANGLES * sizeof(float), cudaHostAllocMapped);
//cudaHostGetDevicePointer(&d_solarAngles, sAnalysis.solarAngles, 0);