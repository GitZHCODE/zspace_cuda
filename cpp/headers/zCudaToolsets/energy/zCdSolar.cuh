// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2022 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : 
//

#ifndef ZSPACE_CD_SOLAR
#define ZSPACE_CD_SOLAR

#pragma once

#include<headers/zCudaToolsets/base/zCdUtilities.cuh>
#include<headers/zCudaToolsets/energy/zTsSolar.h>

using namespace zSpace;

//---- DEVICE VARIABLES

float *d_norms_sunVecs;
float *d_colors;
float *d_cumulativeRadiation;

float* d_hotspotNormals;
float* d_hotspotPositions;
float* d_occlusion;  
zBuildingScore* d_buildingScore;

float* d_buffer_facePositions;
int* d_bufferValues;

int d_solarMemSize;
int d_hotspotMemSize;
int d_occlusionMemSize;


//---- MEMORY METHOD 

ZSPACE_CUDA_EXTERN void cleanSolarDeviceMemory()
{
	cudaFree(d_norms_sunVecs);
	cudaFree(d_colors);
	cudaFree(d_cumulativeRadiation);
}

ZSPACE_CUDA_EXTERN void cleanHotspotDeviceMemory()
{
	cudaFree(d_hotspotNormals);
	cudaFree(d_hotspotPositions);
	cudaFree(d_occlusion);
	cudaFree(d_buildingScore);
}

ZSPACE_CUDA_EXTERN void cleanOcclusionDeviceMemory()
{
	cudaFree(d_buffer_facePositions);
	cudaFree(d_bufferValues);
}

ZSPACE_CUDA_CALLABLE_HOST void setSolarDeviceMemory(int _numNormals, int _numSunVecs)
{
	if (_numNormals + _numSunVecs < d_solarMemSize) return;
	else
	{
		while (d_solarMemSize < _numNormals + _numSunVecs) d_solarMemSize += d_MEMORYMULTIPLIER;

		cleanSolarDeviceMemory();

		checkCudaErrors(cudaMalloc((void**)&d_norms_sunVecs, d_solarMemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_cumulativeRadiation, d_solarMemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_colors, _numNormals * FloatSize));
	}
}

ZSPACE_CUDA_CALLABLE_HOST void setHotspotDeviceMemory(int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals)
{
	if (_numBuildingNormals + _numSunNormals + _numReceiverNormals < d_hotspotMemSize) return;
	else
	{
		while (d_hotspotMemSize < _numBuildingNormals + _numSunNormals + _numReceiverNormals) d_hotspotMemSize += d_MEMORYMULTIPLIER;

		cleanHotspotDeviceMemory();

		checkCudaErrors(cudaMalloc((void**)&d_hotspotNormals, d_hotspotMemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_hotspotPositions, d_hotspotMemSize * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_occlusion, (_numBuildingNormals  / 3) * FloatSize));
		checkCudaErrors(cudaMalloc((void**)&d_buildingScore, (_numBuildingNormals / 3) * sizeof(zBuildingScore)));
	}
}

ZSPACE_CUDA_CALLABLE_HOST void setOcclusionDeviceMemory(int _numBuffers, int _numFacePositions)
{
	if (_numBuffers + _numFacePositions < d_occlusionMemSize) return;
	else
	{
		while (d_occlusionMemSize < _numBuffers + _numFacePositions) d_occlusionMemSize += d_MEMORYMULTIPLIER;
	}

	checkCudaErrors(cudaMalloc((void**)&d_buffer_facePositions, d_occlusionMemSize * FloatSize));
	checkCudaErrors(cudaMalloc((void**)&d_bufferValues, _numBuffers * IntSize));
}

//---- CUDA HOST DEVICE METHODS 

ZSPACE_CUDA_CALLABLE float ofMap(float value, float inputMin, float inputMax, float outputMin, float outputMax)
{
	return ((value - inputMin) / (inputMax - inputMin) * (outputMax - outputMin) + outputMin);
}

ZSPACE_CUDA_CALLABLE zVector getSunPosition(zDate& date, zLocation& location)
{
	float LocalTime = date.tm_hour + (date.tm_min / 60.0);

	double JD = date.toJulian();

	double n = JD - 2451545.0;

	float LDeg = (float)fmod((280.460 + 0.9856474 * n), 360.0);
	float gDeg = (float)fmod((357.528 + 0.9856003 * n), 360.0);

	float LambdaDeg = LDeg + 1.915 * sin(gDeg * DEG_TO_RAD) + 0.01997 * sin(2 * gDeg * DEG_TO_RAD);

	float epsilonDeg = 23.439 - 0.0000004 * n;

	float alphaDeg;
	alphaDeg = atan(cos(epsilonDeg * DEG_TO_RAD) * tan(LambdaDeg * DEG_TO_RAD));
	alphaDeg *= RAD_TO_DEG;
	if (cos(LambdaDeg * DEG_TO_RAD) < 0)	alphaDeg += (4 * (atan(1.0) * RAD_TO_DEG));

	float deltaDeg = asin(sin(epsilonDeg * DEG_TO_RAD) * sin(LambdaDeg * DEG_TO_RAD)) * RAD_TO_DEG;

	zDate dZero(date.tm_year, date.tm_mon, date.tm_mday, 0, 0);
	double JDNull = dZero.toJulian();

	float TNull = ((JDNull - 2451545.0) / 36525);
	float T = LocalTime - location.timeZone;

	float thetaGh = 6.697376 + 2400.05134 * TNull + 1.002738 * T;

	float thetaG = (float)fmod(thetaGh * 15.0, 360.0);
	float theta = thetaG + location.longitude;

	float tauDeg = theta - alphaDeg;

	float denom = (cos(tauDeg * DEG_TO_RAD) * sin(location.latitude * DEG_TO_RAD) - tan(deltaDeg * DEG_TO_RAD) * cos(location.latitude * DEG_TO_RAD));
	float aDeg = atan(sin(tauDeg * DEG_TO_RAD) / denom);
	aDeg *= RAD_TO_DEG;
	if (denom < 0) aDeg = aDeg + 180;
	aDeg += 180; //add 180 to azimuth to compute from the north.

	float hDeg = asin(cos(deltaDeg * DEG_TO_RAD) * cos(tauDeg * DEG_TO_RAD) * cos(location.latitude * DEG_TO_RAD) + sin(deltaDeg * DEG_TO_RAD) * sin(location.latitude * DEG_TO_RAD));
	hDeg *= RAD_TO_DEG;

	float valDeg = hDeg + (10.3 / (hDeg + 5.11));
	float RDeg = 1.02 / (tan(valDeg * DEG_TO_RAD));

	float hRDeg = hDeg + (RDeg / 60);

	return zPoint(cos(aDeg * DEG_TO_RAD) * sin(hRDeg * DEG_TO_RAD), cos(aDeg * DEG_TO_RAD) * cos(hRDeg * DEG_TO_RAD), sin(aDeg * DEG_TO_RAD));
}

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

ZSPACE_CUDA_CALLABLE_DEVICE bool pointInTriangle(zVector& pt, zVector& t0, zVector& t1, zVector& t2)
{
	float Area = 0.5 * (-t1.y * t2.x + t0.y * (-t1.x + t2.x) + t0.x * (t1.y - t2.y) + t1.x * t2.y);

	float s = 1 / (2 * Area) * (t0.y * t2.x - t0.x * t2.y + (t2.y - t0.y) * pt.x + (t0.x - t2.x) * pt.y);
	float t = 1 / (2 * Area) * (t0.x * t1.y - t0.y * t1.x + (t0.y - t1.y) * pt.x + (t1.x - t0.x) * pt.y);

	return ((s >= 0) && (t >= 0) && (s + t <= 1));
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

//---- CUDA KERNEL 

ZSPACE_CUDA_GLOBAL void computeCummulativeRadiation_kernel(float *norms_sunvecs, float* occlusion, float *radiation, float *colors, int numNormals, int numSunvecs, zDomainColor domainColor)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	

	if (i < numNormals/3)
	{
		i *= 3;
		zVector norm(norms_sunvecs[i + 0], norms_sunvecs[i + 1], norms_sunvecs[i + 2]);

		uint sunvecs_offset = numNormals - i;

		float angle = 0;
		float rad = 0;
		int count = 0;
		int radCount = 0;
		for (int o = i; o < i + numSunvecs; o+= 3)  //for (int j =  numNormals; j < numNormals+ numSunvecs; j += 3)
		{
			int j = o + sunvecs_offset;
			

			if (norms_sunvecs[j + 0] != INVALID_VAL && norms_sunvecs[j + 1] != INVALID_VAL && norms_sunvecs[j + 2] != INVALID_VAL)
			{

				zVector sVec(norms_sunvecs[j + 0], norms_sunvecs[j + 1], norms_sunvecs[j + 2]);

				float a = norm.angle(sVec);	
				float weight = ofMap(a, 0.0, 180.0, 1.0, 0.0);		

				if (a > 0.001)
				{
					angle += a;
					count++;
				}
				
				if (weight * radiation[j + 2] > 0.001)
				{
					rad += (weight * radiation[j + 2]);
					radCount++;					
				}				
			}			
		}

		angle /= count;
		rad /= radCount;

		//radiation[i + 0] = angle;
		//radiation[i + 1] = rad;
		radiation[i + 0] = angle * occlusion[i / 3];
		radiation[i + 1] = rad * occlusion[i / 3];
		radiation[i + 2] = INVALID_VAL;

		if (rad < RAD_MIN)
		{
			colors[i + 0] = domainColor.min.h;
			colors[i + 1] = domainColor.min.s;
			colors[i + 2] = domainColor.min.v;
		}
		else if (rad >= RAD_MIN && rad <= RAD_MAX)
		{
			colors[i + 0] = ofMap(rad, RAD_MIN, RAD_MAX, domainColor.min.h, domainColor.max.h);
			colors[i + 1] = ofMap(rad, RAD_MIN, RAD_MAX, domainColor.min.s, domainColor.max.s);
			colors[i + 2] = ofMap(rad, RAD_MIN, RAD_MAX, domainColor.min.v, domainColor.max.v);
		}		
		else
		{
			colors[i + 0] = domainColor.max.h;
			colors[i + 1] = domainColor.max.s;
			colors[i + 2] = domainColor.max.v;
		}


	}
	

}

ZSPACE_CUDA_GLOBAL void computeCummulativeAngle_kernel(float *norms_sunvecs, float *colors, int numNormals, int numSunvecs, zDomainColor domainColor)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < numNormals / 3)
	{
		i *= 3;
		zVector norm(norms_sunvecs[i + 0], norms_sunvecs[i + 1], norms_sunvecs[i + 2]);

		uint sunvecs_offset = numNormals - i;

		float angle = 0;
		int count = 0;
		for (int o = i; o < i + numSunvecs; o += 3)
		{
			int j = o + sunvecs_offset;
			zVector sVec(norms_sunvecs[j + 0], norms_sunvecs[j + 1], norms_sunvecs[j + 2]);

			if (sVec.x != INVALID_VAL && sVec.y != INVALID_VAL && sVec.z != INVALID_VAL)
			{
				float a = norm.angle(sVec);
				if (a > 0.001)
				{
					angle += a;
					count++;
				}				
			}
		}
		angle /= count;



		if (angle > RAD_ANGLE_MAX)
		{
			colors[i + 0] = domainColor.min.h;
			colors[i + 1] = domainColor.min.s;
			colors[i + 2] = domainColor.min.v;
		}
		else if (angle >= RAD_ANGLE_MIN && angle <= RAD_ANGLE_MAX)
		{
			colors[i + 0] = ofMap(angle, RAD_ANGLE_MAX, RAD_ANGLE_MIN, domainColor.min.h, domainColor.max.h);
			colors[i + 1] = ofMap(angle, RAD_ANGLE_MAX, RAD_ANGLE_MIN, domainColor.min.s, domainColor.max.s);
			colors[i + 2] = ofMap(angle, RAD_ANGLE_MAX, RAD_ANGLE_MIN, domainColor.min.v, domainColor.max.v);
		}
		else
		{
			colors[i + 0] = domainColor.max.h;
			colors[i + 1] = domainColor.max.s;
			colors[i + 2] = domainColor.max.v;
    }

		//printf("\n %1.2f %1.2f %1.2f  | %1.2f  | %1.2f %1.2f %1.2f", norm.x, norm.y, norm.z, angle, colors[i + 0], colors[i + 1], colors[i + 2]);



	}


}

ZSPACE_CUDA_GLOBAL void computeOcclusion_kernel(float* buffer_facepositions, int* buffervalues, int numBuffer, int numFacePositions)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < numBuffer / 3)
	{
		int t = i * 3;
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
	}
}

ZSPACE_CUDA_GLOBAL void computeHotspot_kernel(float* normals, float* positions, float* occlusion, zBuildingScore* buildingScore, int numReceiverNormals, int numBuildingNormals, 
	int numSunNormals, float angleTolerance, float disTolerance)
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;

	zVector rVec;
	zVector nVec;
	zVector sVec;

	if (i < numBuildingNormals / 3)
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

		uint receiver_offset = numBuildingNormals + numSunNormals - j;

		float minDis = 100000;

		for (int r = j; r < j + numReceiverNormals; r += 3)
		{
			int k = r + receiver_offset;
			int receiverID = floorf((k - numSunNormals - numBuildingNormals) / 3);

			if(occlusion[buildID] == 0) continue;

			rVec.x = positions[k + 0] - positions[i + 0];
			rVec.y = positions[k + 1] - positions[i + 1];
			rVec.z = positions[k + 2] - positions[i + 2];
			
			if (!isReflectedVector(sVec, nVec, rVec, angleTolerance)) continue;
			if (rVec.length() >= minDis) continue;
			if (rVec.length2() == 0) continue;

			minDis = rVec.length();
		}

		//  -------------------------------- optimize based on buildingScroe array  -----------------------------
		for (int r = j; r < j + numReceiverNormals; r += 3)
		{
			int k = r + receiver_offset;
			int receiverID = floorf((k - numSunNormals - numBuildingNormals) / 3);

			if (occlusion[buildID] == 0) continue;

			rVec.x = positions[k + 0] - positions[i + 0];
			rVec.y = positions[k + 1] - positions[i + 1];
			rVec.z = positions[k + 2] - positions[i + 2];

			if (!isReflectedVector(sVec, nVec, rVec, angleTolerance)) continue;
			if (rVec.length() > minDis * (1 + disTolerance)) continue;
			if (rVec.length2() == 0) continue;

			if (buildingScore[buildID].receiverCount >= MAX_RECEIVER - 1) continue;

			buildingScore[buildID].receiverID[buildingScore[buildID].receiverCount] = receiverID;
			buildingScore[buildID].receiverCount++;
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

ZSPACE_CUDA_EXTERN bool cdpCummulativeRadiation(zTsSolar& solar, float* occlusion, bool EPWRead)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_NORMALS = solar.numSolarNormals();
	int NUM_SUNVECS = solar.numSolarSunVecs();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	zDomainColor dColor = solar.getDomain_Colors();
	setSolarDeviceMemory(NUM_NORMALS, NUM_SUNVECS);

	cudaEventRecord(start);

	// transfer memory to device
	//checkCudaErrors(cudaMemcpy(d_occlusion, solar.getRawOcclusion(), (d_solarMemSize / 3) * FloatSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_norms_sunVecs, solar.getRawNormals_SunVectors(), d_solarMemSize * FloatSize, cudaMemcpyHostToDevice));
	if (EPWRead) checkCudaErrors(cudaMemcpy(d_cumulativeRadiation, solar.getRawCummulativeRadiation(), d_solarMemSize * FloatSize, cudaMemcpyHostToDevice));

	// Launch Kernel
	//printf("\n Launching CDP kernel to compute solar radiation \n ");
	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_solarMemSize / (double)block.x));

	if (EPWRead) computeCummulativeRadiation_kernel << < grid, block >> > (d_norms_sunVecs, occlusion, d_cumulativeRadiation, d_colors, NUM_NORMALS, NUM_SUNVECS, dColor);
	else computeCummulativeAngle_kernel << < grid, block >> > (d_norms_sunVecs, d_colors, NUM_NORMALS, NUM_SUNVECS, dColor);

	checkCudaErrors(cudaGetLastError());

	// transfer memory to host

	checkCudaErrors(cudaMemcpy(solar.getRawColors(), d_colors, NUM_NORMALS * FloatSize, cudaMemcpyDeviceToHost));
	if (EPWRead) checkCudaErrors(cudaMemcpy(solar.getRawCummulativeRadiation(), d_cumulativeRadiation, d_solarMemSize * FloatSize, cudaMemcpyDeviceToHost));

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("\n gpu %1.8f ms \n", milliseconds);
	printf("\n gpu  d_MemSize %i \n", d_solarMemSize);

	return true;
}

ZSPACE_CUDA_EXTERN bool cdpComputeHotspot(zTsSolar& solar, float angleTolerance, float disTolerance)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_RECEIVER_NORMALS = solar.getNumReceiverNormals();
	int NUM_BUILDING_NORMALS = solar.getNumBuildingNormals();
	int NUM_SUN_NORMALS = solar.getNumSunNormals();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);

	setHotspotDeviceMemory(NUM_RECEIVER_NORMALS, NUM_BUILDING_NORMALS, NUM_SUN_NORMALS);

	// ------------------- -------------  transfer memory to device  -------------------------------
	checkCudaErrors(cudaMemcpy(d_hotspotNormals, solar.getHotspotRawNormals(), d_hotspotMemSize * FloatSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_hotspotPositions, solar.getHotspotRawPositions(), d_hotspotMemSize * FloatSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_occlusion, solar.getRawOcclusion(), (NUM_BUILDING_NORMALS / 3) * FloatSize, cudaMemcpyHostToDevice));

	// ------------------- -------------  Launch Kernel	  -------------------------------
	

	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_hotspotMemSize / (double)block.x));

	computeHotspot_kernel << < grid, block >> > (d_hotspotNormals, d_hotspotPositions, d_occlusion, d_buildingScore, NUM_RECEIVER_NORMALS, NUM_BUILDING_NORMALS, NUM_SUN_NORMALS, angleTolerance, disTolerance);

	checkCudaErrors(cudaGetLastError());

	// ------------------- -------------   transfer memory to host	  -------------------------------		
	checkCudaErrors(cudaMemcpy(solar.getRawBuildingScore(), d_buildingScore, (NUM_BUILDING_NORMALS / 3) * sizeof(zBuildingScore), cudaMemcpyDeviceToHost));

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("\n gpu Total %1.8f ms \n", milliseconds);


	return true;
}

ZSPACE_CUDA_EXTERN bool cdpComputeOcclusion(zTsSolar& solar)
{
	int numSMs, numTB;
	cdpGetAttributes(numSMs, numTB);

	// Allocate device memory
	int NUM_BUFFER = solar.numBuffer();
	int NUM_FACEPOSITIONS = solar.numFacePositions();
	

	cudaEvent_t start, stop;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	setOcclusionDeviceMemory(NUM_BUFFER, NUM_FACEPOSITIONS);

	cudaEventRecord(start);

	// transfer memory to device	
	checkCudaErrors(cudaMemcpy(d_buffer_facePositions, solar.getRawBuffer_FacePositions(), d_occlusionMemSize * FloatSize, cudaMemcpyHostToDevice));
	
	// Launch Kernel	
	dim3 block(d_THREADSPERBLOCK);
	dim3 grid((uint)ceil(d_occlusionMemSize / (double)block.x));

	computeOcclusion_kernel << < grid, block >> > (d_buffer_facePositions, d_bufferValues, NUM_BUFFER, NUM_FACEPOSITIONS);
	
	checkCudaErrors(cudaGetLastError());

	// transfer memory to host
	
	checkCudaErrors(cudaMemcpy(solar.getRawBufferValues(), d_bufferValues, NUM_BUFFER * IntSize, cudaMemcpyDeviceToHost));
	
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("\n gpu %1.8f ms \n", milliseconds);	
	printf("\n gpu d_occlusionMemSize %i \n", d_occlusionMemSize);

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