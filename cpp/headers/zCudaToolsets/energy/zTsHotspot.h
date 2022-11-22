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

#ifndef ZSPACE_TS_HOTSPOT
#define ZSPACE_TS_HOTSPOT

#pragma once

#define MAX_SUNVECS_HOUR 366 * 24 * 3
#define MAX_SUNVECS_DAY 24 * 7 * 3
#define INVALID_VAL -10000.0

#define MAX_RECEIVER 1024* 4
//#define MAX_OPTIMIEZD_RECEIVER 1024 * 3


#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zDomain.h>
#include<headers/zCore/utilities/zUtilsCore.h>

namespace zSpace
{
	struct zBuildingScore
	{
		int receiverCount = 0;
		int receiverID[MAX_RECEIVER];		

		//int optimizedReceiverCount = 0;
		//int optimizedReceiverID[MAX_OPTIMIEZD_RECEIVER];
	};

	/*
	struct zBuildingCurrentSunScore
	{
		int receiverCount = 0;
		int receiverID[MAX_SUN_RECEIVER];
	};

	struct zOpBuildingCurrentSunScore
	{
		int receiverCount = 0;
		int receiverID[MAX_OP_SUN_RECEIVER];
	};
	*/


	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsHotspot
	*	\brief A tool set to do hotspot analysis.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_CUDA_TOOLSET zTsHotspot
	{
	private:

		/*!	\brief pointer to container of vertex normals of Building & sunVect & Receiver */
		float* normals;

		/*!	\brief pointer to container of vertex positions of Building & sunVect & Receiver */
		float* positions;

		/*!	\brief number of suns vectors in the container */
		int numSunNormals;

		/*!	\brief number of buildings' normals in the container */
		int numBuildingNormals;

		/*!	\brief number of receivers' normals in the container */
		int numReceiverNormals;

		/*!	\brief pointer to container of scores based on all indexes of receiver & sun  */
		float* receiverScores;

		/*!	\brief pointer to index of sun & building of each receiver */
		//int* scoreIndex;	

		/*!	\brief pointer to container of occlusion */
		float* occlusion;

		/*!	\brief pointer to container of struct buildingScore, for GPU compoutation */
		zBuildingScore* buildingScore;

		/*!	\brief memory size of container of normals */
		int memSize;

	public:

		/*!	\brief compute time of GPU method */
		float computeTime;

		/*!	\brief time of copy memory from GPU to CPU */
		float copyTime;

	public:

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*	\since version 0.0.4
		*/

		ZSPACE_CUDA_CALLABLE_HOST zTsHotspot();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/		

		ZSPACE_CUDA_CALLABLE_HOST ~zTsHotspot();

		//--------------------------
		//---- SET METHODS
		//--------------------------			

		/*! \brief This method sets the normals array
		*
		*	\param		[in]	_receiverNormals					- pointer to container of receivers's normals
		*	\param		[in]	_BuildingNormals					- pointer to container of buildings's normals
		*	\param		[in]	_SunNormals					        - pointer to container of suns's normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setNormals(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals);

		/*! \brief This method sets the positions array
		*
		*	\param		[in]	_receiverPos					    - pointer to container of receivers's positions
		*	\param		[in]	_BuildingPos					    - pointer to container of buildings's positions
		*	\param		[in]	_SunPos					            - pointer to container of suns's positions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setPostion(const float* _receiverPos, const float* _BuildingPos, const float* _SunPos);

		/*! \brief This method sets the array of scores, size of array = numReceiverNormals * numSunNormals / 3
		*
		*	\param		[in]	_receiverScores				    	- pointer to container of scores based on all indexes in the container of receiver & sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setReceiverScores(const float* _receiverScores);	

		/*! \brief This method sets occlusion array , size of array = numBuldingNormals * numSunNormals / 3
		*
		*	\param		[in]	_occlusion					        - pointer to container of occlusion based on all indexes in the container of building & sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setOcclusion(const float* _occlusion);

		/*! \brief This method sets all arrays and numbers that are needed for hotspot computation
		*
		*	\param		[in]	_receiverNormals					- pointer to container of receivers's normals
		*	\param		[in]	_BuildingNormals					- pointer to container of buildings's normals
		*	\param		[in]	_SunNormals					        - pointer to container of suns's normals
		*	\param		[in]	_receiverPos					    - pointer to container of receivers's positions
		*	\param		[in]	_BuildingPos					    - pointer to container of buildings's positions
		*	\param		[in]	_SunPos					            - pointer to container of suns's positions
		*  	\param		[in]	_numReceiverNormals				    - size of container of receivers' normals
		*	\param		[in]	_numBuildingNormals					- size of container of buildings' normals
		*	\param		[in]	_numSunNormals					    - size of container of suns' normals
		*  	\param		[in]	_receiverScores				    	- pointer to container of scores based on all indexes in the container of receiver & sun vectors
		*	\param		[in]	_occlusion					        - pointer to container of occlusion based on all indexes in the container of building & sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspot
		(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals,const float* _receiverPos, const float* _BuildingPos, const float* _SunPos,
			int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals,const float* _receiverScores, const float* _occlusion);

		//--------------------------
		//---- GET METHODS
		//--------------------------	

		/*! \brief This method returns number of sun vectors
		*   \returns int  - number of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int getNumSunNormals();


		/*! \brief This method returns number of buildings' normals
		*   \returns int  - number of buildings' normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int getNumBuildingNormals();

		/*! \brief This method returns number of receivers' normals
		*   \returns int  - number of receivers' normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int getNumReceiverNormals();

		/*! \brief This method returns the pointer to container of normals
		*   \returns float*  - pointer to container of normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawNormals();

		/*! \brief This method returns the pointer to container of positions
		*   \returns float*  - pointer to container of positions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawPositions();

		/*! \brief This method returns the pointer to container of receivers' scores based on the input index in the container of sun vectors
		*  	\param		[in]	sunID       					    - the index of current sun vector
		*   \returns float*  - pointer to container of receivers' scores based on the input index in the container of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getSunIDReceiverScores(int sunID);

		/*! \brief This method returns the pointer to container of receivers' scores based on all indexes in the container of sun vectors
		*   \return float*  - pointer to container of receivers' scores based on all indexes in the container of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getAllReceiverScores();		

		/*! \brief This method computes the highest score of receiver based on the input index in the container of sun vectors
		*  	\param		[in]	sunID       					    - the index of current sun vector in the container
		*   \returns float  - the highest score of receiver based on the input index in the container of sun vectors
		*	\since version 0.0.4 
		*/
		ZSPACE_CUDA_CALLABLE float getSunIDHighScore(int sunID);

		/*! \brief this method returns the highest receiver's score based on all indexes in the container of sun vectors 
		*   \returns float  - the highest receiver's score based on all indexes in the container of sun vectors 
		*	\since version 0.0.4 
		*/
		ZSPACE_CUDA_CALLABLE float getAllHighScore();

		/*! \brief This method computes and returns the sum of all receivers' scores 
		*   \returns float  - the sum of all receivers' scores
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float getAllScoreSum();

		/*! \brief This method returns the pointer to container of occlusion 
		*   \returns float*  - the pointer to container of occlusion 
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawOcclusion();

		/*! \brief This method returns the pointer to container of buildingScore, for GPU computation 
		*   \returns zBuildingScore*  - the pointer to container of buildingScore, for GPU computation
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zBuildingScore* getRawBuildingScore();

		/*! \brief This method cleans the previous computed data of zTsHotspot before another computation 
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void cleanScore();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------	

		/*! \brief This method determines whether given three input vectors satisfy the reflection condition within the tolerance range
		*	\param		[in]	incidenceVector					    - input incidence vector, its direction is from sun to building
		*  	\param		[in]	normalVector					    - input normal vector, its direction is same as vertex normal
		*  	\param		[in]	reflectVector					    - input Reflection vector to be check, its direction is from building to receiver
		*  	\param		[in]	angleTolerance					    - angle tolerance with range of [0.0, 1.0], 0.0 means no tolerance
		*   \returns bool  - true if given three vectors satisfy the reflection condition within the tolerance range, else false.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE bool isReflectedVector(zVector incidenceVector, zVector normalVector, zVector reflectVector, float angleTolerance);

		/*! \brief This method computes the receivers' score in CPU
		*	\param		[in]	angleTolerance					    - angle tolerance with range of [0.0, 0.1], 0.0 means no tolerance
		*  	\param		[in]	disTolerance					    - distance tolerance for optimization, with range of [0.0, 1.0], 0.0 means no tolerance
		*  	\param		[in]	B_occlusion					        - true if enable occlusion condition, else false
		*  	\param		[in]	optimization					    - true if enable optimization method, else false
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeHotspot(float angleTolerance, float disTolerance, bool B_occlusion, bool optimization);

		/*! \brief This method computes receivers' scores based on GPU output data, only compute based on current sun index
		*  	\param		[in]	optimization					    - true if enable optimization method, else false
		*  	\param		[in]	sunID       					    - the index of current sun
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeReceiverScore();
		
		
		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------		

		/*! \brief This method sets memories for containers in zTsHotspot */
		ZSPACE_CUDA_CALLABLE_HOST void setMemory();

	};
}

#if defined(ZSPACE_CUDA_STATIC_LIBRARY)  || defined(ZSPACE_CUDA_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/energy/zTsHotspot.cpp>
#endif

#endif