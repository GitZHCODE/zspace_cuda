

#ifndef ZSPACE_CD_HOTSPOT
#define ZSPACE_CD_HOTSPOT

#pragma once


#include<headers/zCudaToolsets/base/zCdUtilities.cuh>
#include <headers/zCudaToolsets/energy/zTsHotspot.h>

using namespace zSpace;

//---- DEVICE VARIABLES


float* d_normals;
float* d_positions;
float* d_occlusion;  // _numBuildingNormals * _numSunNormals / 3
zBuildingScore* d_buildingScore;  //  numBuildingNormals / 3
int d_MemSize;  // _numReceiverNormals + _numBuildingNormals + _numSunNormals

//---- CUDA HOST DEVICE METHODS 
ZSPACE_CUDA_CALLABLE bool isReflectedVector(zVector incidenceVector, zVector normalVector, zVector reflectVector, float angleTolerance)
{
	float length = (reflectVector.length() > 1) ? reflectVector.length() : 1;

	incidenceVector.normalize();
	normalVector.normalize();
	reflectVector.normalize();

	zVector tureReflectVect = incidenceVector - normalVector * (incidenceVector * normalVector) * 2;

	if ((tureReflectVect * reflectVector) >= (1 - (angleTolerance / length))
		&& (incidenceVector * normalVector) < 0
		&& (reflectVector * normalVector) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}

ZSPACE_CUDA_EXTERN void cleanHotspotDeviceMemory(bool optimization)
{
	// Free memory.
	cudaFree(d_normals);
	cudaFree(d_positions);
	cudaFree(d_occlusion);	
	cudaFree(d_buildingScore);
}

ZSPACE_CUDA_CALLABLE_HOST void setHotspotDeviceMemory(int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals,bool optimization)
{
	if (_numBuildingNormals + _numSunNormals + _numReceiverNormals < d_MemSize) return;
	else
	{
		while (d_MemSize < _numBuildingNormals + _numSunNormals + _numReceiverNormals) d_MemSize += d_MEMORYMULTIPLIER;
	}

	cleanHotspotDeviceMemory(optimization);

	checkCudaErrors(cudaMalloc((void**)&d_normals, d_MemSize * FloatSize));
	checkCudaErrors(cudaMalloc((void**)&d_positions, d_MemSize * FloatSize));
	checkCudaErrors(cudaMalloc((void**)&d_occlusion, (_numBuildingNormals * _numSunNormals / 3) * FloatSize));	
	checkCudaErrors(cudaMalloc((void**)&d_buildingScore, (_numBuildingNormals/3) * sizeof(zBuildingScore)));
}

//---- CUDA KERNEL 

ZSPACE_CUDA_GLOBAL void computeHotspot_kernel
(float* normals, float* positions, float* occlusion,  zBuildingScore* buildingScore,
	int numReceiverNormals, int numBuildingNormals, int numSunNormals, 
	float angleTolerance, float disTolerance, bool B_occlusion, bool optimization)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	
	zVector rVec;
	zVector nVec;
	zVector sVec;

	if (i < numBuildingNormals/3) 
	{
		i *= 3;
		uint sun_offset = numBuildingNormals - i;
		
		nVec.x = normals[i + 0];
		nVec.y = normals[i + 1];
		nVec.z = normals[i + 2];
		int buildID = floorf(i / 3);

		buildingScore[buildID].receiverCount = 0;		

		int j = i + sun_offset;
		int sunID = 0;	
		sVec.x = normals[j + 0];
		sVec.y = normals[j + 1];
		sVec.z = normals[j + 2];

		//__shared__ float share_sth[1000];
		//__syncthreads();

		uint receiver_offset = numBuildingNormals + numSunNormals - j;

		float minDis = 100000;

		for (int r = j; r < j + numReceiverNormals; r += 3)
		{
			int k = r + receiver_offset;
			int receiverID = floorf((k - numSunNormals - numBuildingNormals) / 3);


			rVec.x = positions[k + 0] - positions[i + 0];
			rVec.y = positions[k + 1] - positions[i + 1];
			rVec.z = positions[k + 2] - positions[i + 2];

			if (!isReflectedVector(sVec, nVec, rVec, angleTolerance)) continue;
			if (rVec.length() >= minDis) continue;
			if (rVec.length2() == 0) continue;

			
			//if ((rVec.x * normals[k] + rVec.y * normals[k+1] + rVec.z * normals[k+2] )< 0)continue;
			

			minDis = rVec.length();

			/*
			if (rVec.length2() != 0 && isReflectedVector(sVec, nVec, rVec, angleTolerance))
			{
				if (rVec.length() < minDis)
				{
					minDis = rVec.length();
				}
			}
			*/
		}

		//  -------------------------------- optimize based on buildingScroe array  -----------------------------
		for (int r = j; r < j + numReceiverNormals; r += 3)
		{
			int k = r + receiver_offset;
			int receiverID = floorf((k - numSunNormals - numBuildingNormals) / 3);

			rVec.x = positions[k + 0] - positions[i + 0];
			rVec.y = positions[k + 1] - positions[i + 1];
			rVec.z = positions[k + 2] - positions[i + 2];

			if(!isReflectedVector(sVec, nVec, rVec, angleTolerance)) continue;
			if (rVec.length() > minDis * (1 + disTolerance)) continue;
			if (rVec.length2() == 0) continue;
			
			//if ((rVec.x * normals[k] + rVec.y * normals[k+1] + rVec.z * normals[k+2] )< 0)continue;


			if (buildingScore[buildID].receiverCount >= MAX_RECEIVER - 1) continue;

			buildingScore[buildID].receiverID[buildingScore[buildID].receiverCount] = receiverID;
			buildingScore[buildID].receiverCount++;
			//same subdivison, save time from 1487ms to 685ms; from 49s to 26.7s

			/*
			if (rVec.length2() != 0 && isReflectedVector(sVec, nVec, rVec, angleTolerance) 
				&& rVec.length() <= minDis * (1 + disTolerance))
			{
				buildingScore[buildID].receiverID[buildingScore[buildID].receiverCount] = receiverID;
				buildingScore[buildID].receiverCount++;
			}
			*/
		}

		/*
		if (buildingScore[buildID].receiverCount > 4000 || buildingScore[buildID].optimizedReceiverCount > 1000)
		{
			printf("\n buildID maxCount OP_maxCount = %i  %i %i", buildID, buildingScore[buildID].receiverCount, buildingScore[buildID].optimizedReceiverCount);
		}
		*/
	}
}


//---- launch KERNEL METHODS
ZSPACE_CUDA_EXTERN bool cdpHotspot(zTsHotspot& sHotspot, int sunID, float angleTolerance, float disTolerance, bool B_occlusion, bool optimization)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_RECEIVER_NORMALS = sHotspot.getNumReceiverNormals();
	int NUM_BUILDING_NORMALS = sHotspot.getNumBuildingNormals();
	int NUM_SUN_NORMALS = sHotspot.getNumSunNormals();

	cudaEvent_t start, stop, startKernel,stopKernel, startSet, stopSet;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventCreate(&startKernel);
	cudaEventCreate(&stopKernel);
	cudaEventCreate(&startSet);
	cudaEventCreate(&stopSet);

	cudaEventRecord(start);

	setHotspotDeviceMemory(NUM_RECEIVER_NORMALS, NUM_BUILDING_NORMALS, NUM_SUN_NORMALS, optimization);

	// ------------------- -------------  transfer memory to device  -------------------------------
	checkCudaErrors(cudaMemcpy(d_normals, sHotspot.getRawNormals(), d_MemSize * FloatSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_positions, sHotspot.getRawPositions(), d_MemSize * FloatSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_occlusion, sHotspot.getRawOcclusion(), (NUM_BUILDING_NORMALS * NUM_SUN_NORMALS / 3) * FloatSize, cudaMemcpyHostToDevice));
	
	// ------------------- -------------  Launch Kernel	  -------------------------------
	cudaEventRecord(startKernel);

	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_MemSize / (double)block.x));
	
	computeHotspot_kernel << < grid, block >> > (d_normals, d_positions, d_occlusion, d_buildingScore, NUM_RECEIVER_NORMALS, NUM_BUILDING_NORMALS, NUM_SUN_NORMALS, angleTolerance, disTolerance, B_occlusion, optimization);


	cudaEventRecord(stopKernel);
	cudaEventSynchronize(stopKernel);
	float millisecondsKernel = 0;
	cudaEventElapsedTime(&millisecondsKernel, startKernel, stopKernel);

	checkCudaErrors(cudaGetLastError());

	// ------------------- -------------   transfer memory to host	  -------------------------------		
	checkCudaErrors(cudaMemcpy(sHotspot.getRawBuildingScore(), d_buildingScore, (NUM_BUILDING_NORMALS / 3) * sizeof(zBuildingScore), cudaMemcpyDeviceToHost));
		
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	sHotspot.computeTime = millisecondsKernel;
	sHotspot.copyTime = milliseconds - millisecondsKernel;

	printf("\n gpu Kernel %1.8f ms \n", millisecondsKernel);
	printf("\n gpu Copy %1.8f ms \n", milliseconds - millisecondsKernel);
	printf("\n gpu Total %1.8f ms \n", milliseconds);

	//printf("\n gpu  d_MemSize %i \n", d_MemSize);	
	
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