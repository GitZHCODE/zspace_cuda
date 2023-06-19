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
float* d_bufferValues;

int d_MemSize;

//---- CUDA HOST DEVICE METHODS 


ZSPACE_CUDA_CALLABLE bool pointOnLine(zPoint& pt, zPoint& pA, zPoint& pB, float tolerance = distanceTolerance)
{
	float d = (pt.distanceTo(pA) + pt.distanceTo(pB)) - (pA.distanceTo(pB));
	return (d < tolerance);
}

ZSPACE_CUDA_CALLABLE bool pointInPlanarPolygon(zPoint& pt, zPoint& p0, zPoint& p1, zPoint& p2, zPoint& p3, zVector& planeNormal)
{

	// get component with largest component of normal
	int largeComponent = 0;
	float largeValue = 0;

	if (abs(planeNormal.x) > largeValue)
	{
		largeComponent = 0;
		largeValue = abs(planeNormal.x);
	}

	if (abs(planeNormal.y) > largeValue)
	{
		largeComponent = 1;
		largeValue = abs(planeNormal.y);
	}

	if (abs(planeNormal.z) > largeValue)
	{
		largeComponent = 2;
		largeValue = abs(planeNormal.z);
	}

	// convert points to 2D points by ignorning largest component

	zPoint projectedPt;
	int counter = 0;
	for (int i = 0; i < 3; i++)
	{
		if (i != largeComponent)
		{
			if (counter == 0)
			{
				if (i == 0)  projectedPt.x = pt.x;
				if (i == 1)  projectedPt.x = pt.y;
				if (i == 2)  projectedPt.x = pt.z;


			}
			if (counter == 1)
			{
				if (i == 0) projectedPt.y = pt.x;
				if (i == 1) projectedPt.y = pt.y;
				if (i == 2) projectedPt.y = pt.z;

			}

			counter++;
		}
	}

	zPoint projectedPoints[4];

	for (int j = 0; j < 4; j++)
	{
		zPoint p;
		if (j == 0) p = p0;
		if (j == 1) p = p1;
		if (j == 2) p = p2;
		if (j == 3) p = p3;


		zPoint projP;
		int count = 0;
		for (int i = 0; i < 3; i++)
		{
			if (i != largeComponent)
			{
				if (count == 0)
				{
					if (i == 0)  projP.x = p.x;
					if (i == 1)  projP.x = p.y;
					if (i == 2)  projP.x = p.z;

				}
				if (count == 1)
				{
					if (i == 0)  projP.y = p.x;
					if (i == 1)  projP.y = p.y;
					if (i == 2)  projP.y = p.z;

				}

				count++;
			}


		}

		projP = projP - projectedPt;
		projectedPoints[j] = (projP);

		//cout << endl << points[j];

	}

	projectedPt = zPoint(0, 0, 0);


	// compute winding number
	float windingNum = 0;

	// loop through all edges of the polygon
	//https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf	
	for (int j = 0; j < 4; j++)
	{
		int next = (j + 1) % 4;

		if (projectedPoints[j].y * projectedPoints[next].y < 0)
		{
			float numerator = projectedPoints[j].y * (projectedPoints[next].x - projectedPoints[j].x);
			float denominator = (projectedPoints[j].y - projectedPoints[next].y);
			float r = projectedPoints[j].x + (numerator / denominator);

			if (r > 0)
			{
				if (projectedPoints[j].y < 0) windingNum += 1;
				else windingNum -= 1;
			}
		}
		else if (projectedPoints[j].y == 0 && projectedPoints[j].x > 0)
		{
			if (projectedPoints[next].y > 0) windingNum += 0.5;
			else windingNum -= 0.5;
		}
		else if (projectedPoints[next].y == 0 && projectedPoints[next].x > 0)
		{
			if (projectedPoints[j].y < 0) windingNum += 0.5;
			else windingNum -= 0.5;
		}
	}

	if (windingNum == 0)
	{
		// check if point lies on one of the edges
		for (int j = 0; j < 4; j++)
		{

			int next = (j + 1) % 4;

			float d = (projectedPt.distanceTo(projectedPoints[j]) + projectedPt.distanceTo(projectedPoints[next])) - (projectedPoints[j].distanceTo(projectedPoints[next]));
			
			bool check = ((d < distanceTolerance)) ? true : false; /*pointOnLine(projectedPt, projectedPoints[j], projectedPoints[next])*/;
			if (check)  windingNum = 1;
		}
	}

	if (windingNum != 0)
		return true;
	else
		return false;
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

ZSPACE_CUDA_CALLABLE bool line_lineIntersection(zPoint& a0, zPoint& a1, zPoint& b0, zPoint& b1)
{
	bool out = false;

	zVector u = a1 - a0;
	zVector v = b1 - b0;
	zVector w = a0 - b0;

	double uu = u * u;
	double uv = u * v;
	double vv = v * v;
	double uw = u * w;
	double vw = v * w;

	double denom = ((uu * vv) - (uv * uv));

	if (denom != 0)
	{
		denom = 1 / denom;
		out = true;
	}

	return out;
}

ZSPACE_EXTERN void cleanDeviceMemory()
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
		cleanDeviceMemory();

		checkCudaErrors(cudaMalloc((void **)&d_buffer_facePositions, d_MemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_bufferValues, _numBuffers * FloatSize));
}

//---- CUDA KERNEL 
ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float *buffer_facepositions, float *buffervalues, int numBuffer, int numFacePositions)
{
	uint d_i = blockIdx.x * blockDim.x + threadIdx.x;

	if (d_i < numBuffer && d_i % 3 == 0)
	{
		zPoint a0(buffer_facepositions[d_i + 0], buffer_facepositions[d_i + 1], buffer_facepositions[d_i + 2]);
		zVector dir(10000, 0, 0);
		zPoint a1 = a0 + dir;

		uint faceposition_offset = numBuffer - d_i;

		float dist = 10000000.0;
		float val = -1.0;

		zPoint p0, p1, p2, p3;
		for (int o = d_i; o < d_i + numFacePositions; o += 16)
		{
			int j = o + faceposition_offset;
			float faceID = ((j - numBuffer) / 16);

			p0.x = (buffer_facepositions[j + 0]);  p0.y = (buffer_facepositions[j + 1]);  p0.z = (buffer_facepositions[j + 2]);
			p1.x = (buffer_facepositions[j + 4]);  p1.y = (buffer_facepositions[j + 5]);  p1.z = (buffer_facepositions[j + 6]);
			p2.x = (buffer_facepositions[j + 8]);  p2.y = (buffer_facepositions[j + 9]);  p2.z = (buffer_facepositions[j + 10]);
			p3.x = (buffer_facepositions[j + 12]); p3.y = (buffer_facepositions[j + 13]); p3.z = (buffer_facepositions[j + 14]);

			

			/*check distanceand face winding
			//only compute when current face distance is smaller than stored
			//and the face is in correct winding direction
			*/
			if (buffer_facepositions[j + 3] < dist /*&& !polyWindingXY(p0, p1, p2, p3)*/)
			{
				bool pointInPoly = /*false*/pointInPlanarPolygon(a0, p0, p1, p2, p3, zVector(0, 0, 1));
					//// ASSIGN
							
				if(pointInPoly)
				{
					val = /*-1.0 */faceID;
					dist = buffer_facepositions[j + 3];
				}
			}
		}
		buffervalues[d_i] = val;
		//printf("\n buffervalues %i, %i, %i", i, buffervalues[i], numBuffer);
	}
}

//---- launch KERNEL METHODS

ZSPACE_EXTERN bool cdComputeOcclusion(zTsSolarOcclusionGPU &sOcclusion)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_BUFFER = sOcclusion.numBuffer();
	int NUM_FACEPOSITIONS = sOcclusion.numFacePositions();

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
	computeOcclusion_kernel << < grid, block >> > (d_buffer_facePositions, d_bufferValues, NUM_BUFFER, NUM_FACEPOSITIONS);
	cudaEventRecord(stop_kernel);
	cudaEventSynchronize(stop_kernel);

	checkCudaErrors(cudaGetLastError());

	// transfer memory to host
	cudaEventRecord(start_toHost);
	checkCudaErrors(cudaMemcpy(sOcclusion.getRawBufferValues(), d_bufferValues, NUM_BUFFER * FloatSize, cudaMemcpyDeviceToHost));
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
	printf("\n time_kernel %1.8f ms \n", time_kernel);
	printf("\n time_toHost %1.8f ms \n", time_toHost);

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