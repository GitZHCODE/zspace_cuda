// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2022 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Ling Mao <Ling.Mao@zaha-hadid.com>, Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//

#include <headers/zCudaToolsets/energy/zTsHotspot.h>


namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsHotspot::zTsHotspot() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsHotspot::~zTsHotspot() {}

	//---- SET METHODS

	ZSPACE_INLINE void zTsHotspot::setNormals
	(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals)
	{
		std::copy(_BuildingNormals, _BuildingNormals + numBuildingNormals, normals);
		std::copy(_SunNormals, _SunNormals + numSunNormals, normals + numBuildingNormals);
		std::copy(_receiverNormals, _receiverNormals + numReceiverNormals, normals + numSunNormals + numBuildingNormals);
	}

	ZSPACE_INLINE void zTsHotspot::setPostion
	(const float* _receiverPos, const float* _BuildingPos, const float* _SunPos)
	{	
		std::copy(_BuildingPos, _BuildingPos + numBuildingNormals, positions);
		std::copy(_SunPos, _SunPos + numSunNormals, positions + numBuildingNormals);
		std::copy(_receiverPos, _receiverPos + numReceiverNormals, positions + numSunNormals + numBuildingNormals);
	}

	ZSPACE_INLINE void zTsHotspot::setReceiverScores(const float* _receiverScores)
	{
		std::copy(_receiverScores, _receiverScores + (numReceiverNormals * numSunNormals / 3), receiverScores);
	}

	ZSPACE_INLINE void zTsHotspot::setOcclusion(const float* _occlusion)
	{
		std::copy(_occlusion, _occlusion + (numBuildingNormals * numSunNormals / 3), occlusion);
	}

	ZSPACE_INLINE void zTsHotspot::setHotspot
	(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals,
		const float* _receiverPos, const float* _BuildingPos, const float* _SunPos,
		int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals, 
		const float* _receiverScores, const float* _occlusion)
	{		
		numBuildingNormals = _numBuildingNormals;
		numSunNormals = _numSunNormals;
		numReceiverNormals = _numReceiverNormals;
		computeTime = 0;
		copyTime = 0;

		setMemory();

		setNormals(_receiverNormals, _BuildingNormals, _SunNormals);
		setPostion(_receiverPos, _BuildingPos, _SunPos);
		setReceiverScores(_receiverScores);
		setOcclusion(_occlusion);
	}

	//---- GET METHODS

	ZSPACE_INLINE int zTsHotspot::getNumSunNormals()
	{
		return numSunNormals;
	}

	ZSPACE_INLINE int zTsHotspot::getNumBuildingNormals()
	{
		return numBuildingNormals;
	}

	ZSPACE_INLINE int zTsHotspot::getNumReceiverNormals()
	{
		return numReceiverNormals;
	}

	ZSPACE_INLINE float* zTsHotspot::getRawNormals()
	{
		return normals;
	}

	ZSPACE_INLINE float* zTsHotspot::getRawPositions()
	{
		return positions;
	}

	ZSPACE_INLINE float* zTsHotspot::getSunIDReceiverScores(int sunID)
	{
		float* curReceiverScores = new float[numReceiverNormals];
		std::copy(receiverScores + (sunID * numReceiverNormals), receiverScores+ ((sunID + 1) * numReceiverNormals), curReceiverScores);

		return curReceiverScores;
	}

	ZSPACE_INLINE float* zTsHotspot::getAllReceiverScores()
	{
		return receiverScores;
	}

	ZSPACE_INLINE float zTsHotspot::getSunIDHighScore(int sunID)
	{
		float highestScore = 0;

		for (int i = sunID * numReceiverNormals; i < (sunID + 1) * numReceiverNormals; i += 3)
		{
			if (receiverScores[i] > highestScore)
			{
				highestScore = receiverScores[i];
			}
		}

		return highestScore;
	}

	ZSPACE_INLINE float zTsHotspot::getAllHighScore()
	{
		float highestScore = 0;

		for (int i = 0; i < (numBuildingNormals * numSunNormals / 3); i += 3)
		{
			if (receiverScores[i] > highestScore)
			{
				highestScore = receiverScores[i];
			}
		}

		return highestScore;
	}

	ZSPACE_INLINE float zTsHotspot::getAllScoreSum()
	{
		float sum = 0;
		for (int i = 0; i < numSunNormals; i += 3)
		{
			int sunID = floor(i / 3);
			for (int j = 0; j < numReceiverNormals; j += 3)
			{
				sum += receiverScores[sunID * numReceiverNormals + j];
			}
		}
		return sum;
	}

	ZSPACE_INLINE float* zTsHotspot::getRawOcclusion()
	{
		return occlusion;
	}

	ZSPACE_INLINE zBuildingScore* zTsHotspot::getRawBuildingScore()
	{
		return buildingScore;
	}

	ZSPACE_INLINE void zTsHotspot::cleanScore()
	{
		for (int j = 0; j < numReceiverNormals; j += 3)
		{
			receiverScores[j] = 0;
		}

		for (int i = 0; i < numBuildingNormals; i += 3)
		{
			int buildID = floor(i / 3);
			buildingScore[buildID].receiverCount = 0;
		}

	}

	//---- COMPUTE METHODS

	ZSPACE_INLINE bool zTsHotspot::isReflectedVector(zVector incidenceVector, zVector normalVector, zVector reflectVector, float angleTolerance)
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

	ZSPACE_INLINE void zTsHotspot::computeHotspot(float angleTolerance, float disTolerance, bool B_occlusion, bool optimization)
	{
		
		zVector rVec;
		zVector nVec;
		zVector sVec;
		int maxCount = 0;
		int maxCount2 = 0;

		for (int i = 0; i < numBuildingNormals; i += 3)
		{			
			nVec.x = normals[i + 0];
			nVec.y = normals[i + 1];
			nVec.z = normals[i + 2];

			int buildID = floor(i/3);
			int count = 0;
			int count2 = 0;

			for (int j = numBuildingNormals; j < numSunNormals + numBuildingNormals; j += 3)
			{
				int sunID = floor((j - numBuildingNormals) / 3);				

				if (occlusion[sunID * numBuildingNormals + i] == 0 && B_occlusion) continue;

				sVec.x = normals[j + 0];
				sVec.y = normals[j + 1];
				sVec.z = normals[j + 2];

				float minDis = 100000;
				int goodReceiverID[10000];
				int goodReceiverCount = 0;

				for (int k = numSunNormals + numBuildingNormals; k < numSunNormals + numBuildingNormals + numReceiverNormals; k += 3)
				{
					rVec.x = positions[k + 0] - positions[i + 0];
					rVec.y = positions[k + 1] - positions[i + 1];
					rVec.z = positions[k + 2] - positions[i + 2];

					//printf("%i %i %i \n", rVec.x, rVec.y, rVec.z);					

					if (rVec.length2() != 0 && isReflectedVector(sVec, nVec, rVec, angleTolerance))
					{
						goodReceiverID[goodReceiverCount] = k;
						receiverScores[numReceiverNormals * sunID + k - numSunNormals - numBuildingNormals] ++;
						goodReceiverCount++;
						//count++; 
						//count2++;

						if (rVec.length() < minDis)
						{
							minDis = rVec.length();
						}
					}
				}
				//printf("\n %i  %i  %1.8f ", buildID, sunID, minDis);

				if (optimization)
				{
					for (int r = 0; r < goodReceiverCount; r += 1)
					{
						int k = goodReceiverID[r];

						rVec.x = positions[k + 0] - positions[i + 0];
						rVec.y = positions[k + 1] - positions[i + 1];
						rVec.z = positions[k + 2] - positions[i + 2];

						if (rVec.length() > minDis *(1 + disTolerance) )
						{
							if (receiverScores[numReceiverNormals * sunID + k - numSunNormals - numBuildingNormals] <= 0)
							{
								cout << "ERROR negetive score, receiverID:  " <<( k - numSunNormals - numBuildingNormals) / 3 <<endl;
							}
							else
							{ 								
								receiverScores[numReceiverNormals * sunID + k - numSunNormals - numBuildingNormals] -= 1;
								//count2--;
							}
						}
					}
				}
			}

			//if (count > maxCount)maxCount = count;
			//if (count2 > maxCount2)maxCount2 = count2;
		}

		//printf("\n maxCount =  %i %i ", maxCount, maxCount2);
	}

	ZSPACE_INLINE void zTsHotspot::computeReceiverScore()
	{
		int maxCount = 0;

		for (int i = 0; i < numBuildingNormals; i += 3)
		{
			int buildID = floor(i / 3);

			for (int r = 0; r < buildingScore[buildID].receiverCount; r++)
			{
				int receiverID = buildingScore[buildID].receiverID[r];
				receiverScores[receiverID * 3]++;
			}

			if (maxCount < buildingScore[buildID].receiverCount)
			{
				maxCount = buildingScore[buildID].receiverCount;
			}
		}

		printf("\n maxCount maxOpCount = %i", maxCount);

		/*for (int i = 0; i < numReceiverNormals; i += 3)
		{
			if (receiverScores[i] > 0)
			{
				printf("\n receiverScores = %i", i);
			}
		}*/
	}

	
	//---- PROTECTED METHODS

	ZSPACE_INLINE void zTsHotspot::setMemory()
	{
		if ((numBuildingNormals + numSunNormals + numReceiverNormals) < memSize) return;
		else
		{
			while (memSize < numBuildingNormals + numSunNormals + numReceiverNormals) memSize += d_MEMORYMULTIPLIER;
		}

		normals = new float[memSize];
		positions = new float[memSize];
		occlusion = new float[numBuildingNormals * numSunNormals / 3];
		receiverScores = new float[numReceiverNormals * numSunNormals / 3];
		buildingScore = new zBuildingScore[numBuildingNormals / 3];
	}
}
