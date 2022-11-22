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

#ifndef ZSPACE_CD_SOLAR_ANALYSIS
#define ZSPACE_CD_SOLAR_ANALYSIS

#pragma once




#include<headers/zCudaToolsets/base/zCdUtilities.cuh>
#include <headers/zCudaToolsets/energy/zTsSolarOcclusionGPU.h>


using namespace zSpace;

//---- DEVICE VARIABLES

float *d_buffer_facePositions;
int* d_bufferValues;

int d_MemSize;
int d_SharedMemSize;

//user define
const int MAX_SHAREDMEMSIZE = 48000;
//const int MAX_SHAREDMEMSIZE = 64;

//---- CUDA HOST DEVICE METHODS 

ZSPACE_CUDA_CALLABLE_DEVICE bool pointInTriangle(zVector& pt, zVector& t0, zVector& t1, zVector& t2)
{
	float Area = 0.5*(-t1.y * t2.x + t0.y * (-t1.x + t2.x) + t0.x * (t1.y - t2.y) + t1.x * t2.y);

	float s = 1 / (2 * Area) * (t0.y * t2.x - t0.x * t2.y + (t2.y - t0.y) * pt.x + (t0.x - t2.x) * pt.y);
	float t = 1 / (2 * Area) * (t0.x * t1.y - t0.y * t1.x + (t0.y - t1.y) * pt.x + (t1.x - t0.x) * pt.y);

	return ((s >= 0) && (t >= 0) && (s+t<=1));


	//// Compute vectors        
	//zVector v0 = t2 - t0;
	//zVector	v1 = t1 - t0;
	//zVector	v2 = pt - t0;

	//// Compute dot products
	//double	dot00 = v0 * v0;
	//double	dot01 = v0 * v1;
	//double	dot02 = v0 * v2;
	//double	dot11 = v1 * v1;
	//double	dot12 = v1 * v2;

	//// Compute barycentric coordinates
	//double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	//double 	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	//double	v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	//if (abs(u) < 0.001) u = 0;
	//if (abs(v) < 0.001) v = 0;

	//// round factor to precision 3 
	//double factor = powf(10, 3);

	//u = std::round(u * factor) / factor;
	//v = std::round(v * factor) / factor;

	//// Check if point is in triangle
	//return ((u >= 0) && (v >= 0) && (u + v <= 1));
}

ZSPACE_CUDA_CALLABLE_DEVICE bool polyWindingXY(zPoint& p0, zPoint& p1, zPoint& p2, zPoint& p3)
{
	float winding = 0.0f;

	for (int i = 0; i < 4; i++)
	{
		zPoint previous;
		zPoint current;

			if (i == 0) { previous = p3; current = p0; }
			if (i == 1) { previous = p0; current = p1; }
			if (i == 2) { previous = p1; current = p2; }
			if (i == 3) { previous = p2; current = p3; }

		winding += (previous.x * current.y - current.x * previous.y);
	}
	if (winding > 0) return true;
	else return false;
}

ZSPACE_CUDA_CALLABLE_DEVICE bool pointInBounds(zPoint& pt, zPoint& p0, zPoint& p1, zPoint& p2, zPoint& p3)
{


	zVector minBB = zVector(10000, 10000, 10000);
	zVector maxBB = zVector(-10000, -10000, -10000);

	zPoint inPoint;
	for (int i = 0; i < 4; i++)
	{
		if (i == 0) inPoint = p0;
		if (i == 1) inPoint = p1;
		if (i == 2) inPoint = p2;
		if (i == 3) inPoint = p3;

		if (inPoint.x < minBB.x) minBB.x = inPoint.x;
		if (inPoint.y < minBB.y) minBB.y = inPoint.y;
		if (inPoint.z < minBB.z) minBB.z = inPoint.z;

		if (inPoint.x > maxBB.x) maxBB.x = inPoint.x;
		if (inPoint.y > maxBB.y) maxBB.y = inPoint.y;
		if (inPoint.z > maxBB.z) maxBB.z = inPoint.z;
	}

	if (pt.x > minBB.x && pt.x<maxBB.x && pt.y>minBB.y && pt.y < maxBB.y) return true;
	else return false;
}

ZSPACE_CUDA_EXTERN void cleanDeviceMemory()
{
	// Free memory.
	cudaFree(d_buffer_facePositions);
	cudaFree(d_bufferValues);
}

ZSPACE_CUDA_CALLABLE_HOST void setDeviceMemory(int _numBuffers, int _numFacePositions)
{
	if (_numBuffers + _numFacePositions < d_MemSize) return;
	else
	{
		while (d_MemSize < _numBuffers + _numFacePositions) d_MemSize += d_MEMORYMULTIPLIER;
	}

	if (_numFacePositions < d_SharedMemSize) return;
	else
	{
		while (d_SharedMemSize < _numFacePositions*sizeof(int)) d_SharedMemSize += MAX_SHAREDMEMSIZE;
	}
		cleanDeviceMemory();

		checkCudaErrors(cudaMalloc((void **)&d_buffer_facePositions, d_MemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_bufferValues, _numBuffers * IntSize));
}

//---- CUDA KERNEL 


//ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float* buffer_facepositions, int* buffervalues, int numBuffer, int numFacePositions, int maxNumFace)
//{
//	uint i = blockIdx.x * blockDim.x + threadIdx.x;
//
//
//	//if (i < numBuffer && i % 3 == 0)
//	if (i < numBuffer / 3)
//		//for (uint i = blockIdx.x * blockDim.x + threadIdx.x;
//		//	i < numBuffer / 3;
//		//	i += blockDim.x * gridDim.x)
//	{
//		int t = i * 3;
//		//parallel buffer pixel for computation
//		zPoint bufferPos(buffer_facepositions[t + 0], buffer_facepositions[t + 1], buffer_facepositions[t + 2]);
//
//		uint faceposition_offset = numBuffer - t;
//
//		int numFace = numFacePositions / 16;
//		int currentNumF;
//		int iteration = floorf(numFace / maxNumFace) + 1;
//
//		float dist = 10000000;
//		int val = -1;
//
//		//assign maximum number of faces for 48kb shared memories in one iteration
//		for (int m = 0; m < iteration; m++)
//		{
//			if (m < iteration - 1) currentNumF = maxNumFace;
//			else currentNumF = numFace - maxNumFace * m;
//			
//			//assign shared memories of maximum number of faces
//			extern __shared__ float shared[];
//			int counter = 0;
//			while (counter < currentNumF * 16)
//			{
//				shared[counter] = buffer_facepositions[numBuffer + m * maxNumFace * 16 + counter];
//				counter++;
//			}
//			printf("\n counter %i", counter);
//
//			//synchronise all threads finished copy
//			__syncthreads();
//
//			//traverse all faces in one iteration
//			//for (int o = 0; o < currentNumF; o++)
//			int o = 0;
//			do
//			{
//				int j = o * 16;
//				int faceID = maxNumFace * m + o;
//
//				zPoint p0(shared[j + 0], shared[j + 1], shared[j + 2]);
//				zPoint p1(shared[j + 4], shared[j + 5], shared[j + 6]);
//				zPoint p2(shared[j + 8], shared[j + 9], shared[j + 10]);
//				zPoint p3(shared[j + 12], shared[j + 13], shared[j + 14]);
//
//				/*check distanceand face winding
//				only compute when current face distance is smaller than stored
//				and the face is in correct winding direction
//				*/
//				//if (pointInBounds(bufferPos, p0, p1, p2, p3))
//				{
//					if (shared[j + 3] < dist && !polyWindingXY(p0, p1, p2, p3))
//					{
//						if (pointInTriangle(bufferPos, p0, p1, p2) || pointInTriangle(bufferPos, p0, p2, p3))
//						{
//							val = faceID;
//							dist = shared[j + 3];
//						}
//					}
//					faceID++;
//					o++;
//				}
//			}
//			while (o < currentNumF);
//		}
//		//only assign buffervalue once after all faces computed
//		buffervalues[i] = val;
//		//printf("\n buffervalues %i, %i, %i", i, buffervalues[i], numBuffer);
//	}
//}

//
//ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float* buffer_facepositions, int* buffervalues, int numBuffer, int numFacePositions, int maxNumFace)
//{
//	uint i = blockIdx.x * blockDim.x + threadIdx.x;
//
//	//if (i < numBuffer && i % 3 == 0)
//	if (i < numBuffer / 3)
//	{
//		int t = i * 3;
//		zPoint bufferPos(buffer_facepositions[t + 0], buffer_facepositions[t + 1], buffer_facepositions[t + 2]);
//
//		uint faceposition_offset = numBuffer - t;
//
//		float dist = 10000000;
//		int val = -1;
//
//		for (int o = t; o < t + numFacePositions; o += 16)
//		{
//			int j = o + faceposition_offset;
//			//int faceID = floorf((j - numBuffer) / 16);
//			int faceID = (j - numBuffer) / 16;
//
//			extern __shared__ float shared[];
//			for (int i = 0; i < 16; i++)
//			{
//				shared[i] = buffer_facepositions[j + i];
//			}
//			__syncthreads();
//			zPoint p0(shared[0], shared[1], shared[2]);
//			zPoint p1(shared[4], shared[5], shared[6]);
//			zPoint p2(shared[8], shared[9], shared[10]);
//			zPoint p3(shared[12], shared[13], shared[14]);
//
//			/*check distanceand face winding
//			only compute when current face distance is smaller than stored
//			and the face is in correct winding direction
//			*/
//			//if (buffer_facepositions[j + 3] < dist && !polyWindingXY(p0, p1, p2, p3))
//			if (pointInBounds(bufferPos, p0, p1, p2, p3))
//			{
//				if (shared[3] < dist && !polyWindingXY(p0, p1, p2, p3))
//				{
//					if (pointInTriangle(bufferPos, p0, p1, p2) || pointInTriangle(bufferPos, p0, p2, p3))
//					{
//						val = faceID;
//						dist = shared[3];
//					}
//				}
//			}
//		}
//		buffervalues[i] = val;
//		//printf("\n buffervalues %i, %i, %i", i, buffervalues[i], numBuffer);
//	}
//}

//
//ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float* buffer_facepositions, int* buffervalues, int numBuffer, int numFacePositions, int maxNumFace)
//{
//	uint i = blockIdx.x * blockDim.x + threadIdx.x;
//
//	__shared__ float shared[16][16];
//
//	if (i < numBuffer / 3)
//	{
//		int t = i * 3;
//		zPoint bufferPos(buffer_facepositions[t + 0], buffer_facepositions[t + 1], buffer_facepositions[t + 2]);
//
//		uint faceposition_offset = numBuffer - t;
//
//		float dist = 10000000;
//		int val = -1;
//
//		for (int o = t; o < t + numFacePositions; o += 16)
//		{
//			int j = o + faceposition_offset;
//			//int faceID = floorf((j - numBuffer) / 16);
//			int faceID = (j - numBuffer) / 16;
//
//			zPoint p0(buffer_facepositions[j + 0], buffer_facepositions[j + 1], buffer_facepositions[j + 2]);
//			zPoint p1(buffer_facepositions[j + 4], buffer_facepositions[j + 5], buffer_facepositions[j + 6]);
//			zPoint p2(buffer_facepositions[j + 8], buffer_facepositions[j + 9], buffer_facepositions[j + 10]);
//			zPoint p3(buffer_facepositions[j + 12], buffer_facepositions[j + 13], buffer_facepositions[j + 14]);
//
//			/*check distanceand face winding
//			only compute when current face distance is smaller than stored
//			and the face is in correct winding direction
//			*/
//			if (pointInBounds(bufferPos, p0, p1, p2, p3))
//			{
//				if (buffer_facepositions[j + 3] < dist && !polyWindingXY(p0, p1, p2, p3))
//				{
//					if (pointInTriangle(bufferPos, p0, p1, p2) || pointInTriangle(bufferPos, p0, p2, p3))
//					{
//						val = faceID;
//						dist = buffer_facepositions[j + 3];
//					}
//				}
//			}
//		}
//		buffervalues[i] = val;
//		//printf("\n buffervalues %i, %i, %i", i, buffervalues[i], numBuffer);
//	}
//}


ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float* buffer_facepositions, int* buffervalues, int numBuffer, int numFacePositions, int maxNumFace)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;

	//if (i < numBuffer && i % 3 == 0)
		if (i < numBuffer/3)
	{
		int t =  i * 3;
		zPoint bufferPos(buffer_facepositions[t + 0], buffer_facepositions[t + 1], buffer_facepositions[t + 2]);

		uint faceposition_offset = numBuffer - t;

		float dist = 10000000;
		int val = -1;

		for (int o = t; o < t + numFacePositions; o += 16)
		{
			int j = o + faceposition_offset;
			//int faceID = floorf((j - numBuffer) / 16);
			int faceID = (j - numBuffer) / 16;

			zPoint p0(buffer_facepositions[j + 0], buffer_facepositions[j + 1], buffer_facepositions[j + 2]);
			zPoint p1(buffer_facepositions[j + 4], buffer_facepositions[j + 5], buffer_facepositions[j + 6]);
			zPoint p2(buffer_facepositions[j + 8], buffer_facepositions[j + 9], buffer_facepositions[j + 10]);
			zPoint p3(buffer_facepositions[j + 12], buffer_facepositions[j + 13], buffer_facepositions[j + 14]);

			/*check distanceand face winding
			only compute when current face distance is smaller than stored
			and the face is in correct winding direction
			*/
			if (pointInBounds(bufferPos, p0, p1, p2, p3))
			{
				if (buffer_facepositions[j + 3] < dist && !polyWindingXY(p0, p1, p2, p3))
				{
					if (pointInTriangle(bufferPos, p0, p1, p2) || pointInTriangle(bufferPos, p0, p2, p3))
					{
						val = faceID;
						dist = buffer_facepositions[j + 3];
					}
				}
			}
		}
		buffervalues[i] = val;
		//printf("\n buffervalues %i, %i, %i", i, buffervalues[i], numBuffer);
	}
}



//---- launch KERNEL METHODS

ZSPACE_CUDA_EXTERN bool cdComputeOcclusion(zTsSolarOcclusionGPU &sOcclusion)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_BUFFER = sOcclusion.numBuffer();
	int NUM_FACEPOSITIONS = sOcclusion.numFacePositions();
	int maxNumFace = MAX_SHAREDMEMSIZE / (sizeof(float) * 16);

	cudaEvent_t
		start, stop,
		start_toDevice, stop_toDevice,
		start_kernel, stop_kernel,
		start_toHost, stop_toHost;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventCreate(&start_toDevice);
	cudaEventCreate(&stop_toDevice);
	cudaEventCreate(&start_kernel);
	cudaEventCreate(&stop_kernel);
	cudaEventCreate(&start_toHost);
	cudaEventCreate(&stop_toHost);
	
	setDeviceMemory(NUM_BUFFER, NUM_FACEPOSITIONS);
  
	cudaEventRecord(start);

	// transfer memory to device
	cudaEventRecord(start_toDevice);
	checkCudaErrors(cudaMemcpy(d_buffer_facePositions, sOcclusion.getRawBuffer_FacePositions(), d_MemSize * FloatSize, cudaMemcpyHostToDevice));
	cudaEventRecord(stop_toDevice);
	cudaEventSynchronize(stop_toDevice);

	// Launch Kernel
	//printf("\n Launching CDP kernel to compute solar radiation \n ");
	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_MemSize / (double)block.x));	
	
	cudaEventRecord(start_kernel);
	//computeOcclusion_kernel << < grid, block >> > (d_buffer_facePositions, d_bufferValues, NUM_BUFFER, NUM_FACEPOSITIONS);
	computeOcclusion_kernel << < grid, block, MAX_SHAREDMEMSIZE >> > (d_buffer_facePositions, d_bufferValues, NUM_BUFFER, NUM_FACEPOSITIONS, maxNumFace);
	cudaEventRecord(stop_kernel);
	cudaEventSynchronize(stop_kernel);

	checkCudaErrors(cudaGetLastError());

	// transfer memory to host
	cudaEventRecord(start_toHost);
	checkCudaErrors(cudaMemcpy(sOcclusion.getRawBufferValues(), d_bufferValues, NUM_BUFFER * IntSize, cudaMemcpyDeviceToHost));
	cudaEventRecord(stop_toHost);
	cudaEventSynchronize(stop_toHost);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	float milliseconds = 0;
	float time_toDevice = 0;
	float time_kernel = 0;
	float time_toHost = 0;
	cudaEventElapsedTime(&time_toDevice, start_toDevice, stop_toDevice);
	cudaEventElapsedTime(&time_kernel, start_kernel, stop_kernel);
	cudaEventElapsedTime(&time_toHost, start_toHost, stop_toHost);

	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("\n gpu %1.8f ms \n", milliseconds);
	printf("\n time_toDevice %1.8f ms \n", time_toDevice);
	printf("time_kernel %1.8f ms \n", time_kernel);
	printf("time_toHost %1.8f ms \n", time_toHost);

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