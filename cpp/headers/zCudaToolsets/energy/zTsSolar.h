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

#ifndef ZSPACE_TS_SOLAR
#define ZSPACE_TS_SOLAR

#pragma once

#define MAX_SUNVECS_HOUR 366*24 * 3
#define MAX_SUNVECS_DAY 24 * 7 * 3
#define MAX_SUNVECS_DATE 3
#define INVALID_VAL -10000.0

#define COMPASS_SUBD 12*2*3

#define RAD_ANGLE_MIN 30.0
#define RAD_ANGLE_MAX 120.0

#define RAD_MIN 300.0
#define RAD_MAX 1000.0

#define MAX_RECEIVER 1024* 4

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
	};

	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsSolar
	*	\brief A tool set to do solar analysis, occlusion and hotspot.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_CUDA_TOOLSET zTsSolar
	{
	private:

		// ------ solar -------
		/*!	\brief pointer to container of building's normals*/
		float* solarNormals;
		
		/*!	\brief number of building's normals in the container*/
		int numSolarNorms;

		/*!	\brief pointer to container of sun vectors for each hour of a year*/
		float* sunVecs_hour;

		/*!	\brief pointer to container of sun vectors for each hour of 7 special days*/
		float* sunVecs_days;

		/*!	\brief pointer to container of sun vectors of a date*/
		float* sunVecs_date;

		/*!	\brief pointer to container of compass display vectors*/
		float* compassPts;

		/*!	\brief pointer to container of epw data*/
		float* epwData_radiation;
		
		/*!	\brief number of epw data points in the container*/
		int numEPWData;

		/*!	\brief date domain*/
		zDomainDate dDate;

		/*!	\brief location information*/
		zLocation location;

		/*!	\brief color domain*/
		zDomainColor dColor;

		/*!	\brief pointer to container of cummulative raditation*/
		float* cummulativeRadiation;

		/*!	\brief pointer to container of vertex's colors based on radiation*/
		float* colors;

		/*!	\brief pointer to container of normals of Building & sunVect, for radiation calculation*/
		float* norm_sunvecs;

		// ------ hotspot -------
		/*!	\brief pointer to container of vertex normals of Building & sunVect & Receiver, for hotspot calculation*/
		float* hotspotNormals;

		/*!	\brief pointer to container of vertex positions of Building & sunVect & Receiver, for hotspot calculation*/
		float* hotspotPositions;
		
		/*!	\brief number of sun vectors in the container */
		int numSunNormals;

		/*!	\brief number of buildings' normals in the container */
		int numBuildingNormals;

		/*!	\brief number of receivers' normals in the container */
		int numReceiverNormals;

		/*!	\brief pointer to container of scores based on all indexes of receivers & suns*/
		float* receiverScores;

		/*!	\brief pointer to container of occlusion*/
		float* occlusion;

		/*!	\brief pointer to container of hotspot score of each building's vertex, for hotspot calculation*/
		zBuildingScore* buildingScore;

		// ------ occlusion -------
		/*!	\brief pointer to container of buffer positions	*/
		float* bufferPositions;

		/*!	\brief pointer to container of buffer values*/
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

		// ------ memory -------
		/*!	\brief memory size of container of solar hotspot */
		int solarMemSize;

		/*!	\brief memory size of container of hotspot  */
		int hotspotMemSize;

		/*!	\brief memory size of container of occlusion  */
		int occlusionMemSize;


	public:

		zUtilsCore coreUtils;

		//float computeTime;

		//float copyTime;
		
		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------
		/*! \brief Default constructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST zTsSolar();
		

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------
		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST ~zTsSolar();

		//-------------------------------------------------------------------------------
		//---- SET METHODS
		
		// ------ solar -------
		/*! \brief This method sets container of normals
		*	\param		[in]	_normals					- pointer to container of normals
		*	\param		[in]	_numNormals					- number of face normals in the container
		*	\param		[in]	EPWread					    - true if reads EPW as input, else false
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setSolarNormals(const float* _normals, int _numNormals, bool EPWread);
		
		/*! \brief This method sets EPW data based on the input path, returns true if read EPW file successfully, else false
		*	\param		[in]	path					    - path of EPW file
		*   \returns bool  - true if read EPW file successfully, else false
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST bool setEPWData(string path);

		/*! \brief This method sets the domain of dates
		*	\param		[in]	_dDate					    - domain of dates
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Dates(zDomainDate& _dDate);
		
		/*! \brief This method sets the domain of colors
		*	\param		[in]	_dColor					    - domain of colors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Colors(zDomainColor& _dColor);

		/*! \brief This method sets the geographical location of site, including timeZone, longitude and latitude
		*	\param		[in]	_location					- geographical location, including timeZone, longitude and latitude
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setLocation(zLocation& _location);

		// ------ hotspot -------

		/*! \brief This method sets the numbers of receivers, buildings and suns
		*
		*	\param		[in]	_numReceiverNormals					- number of receivers's normals
		*	\param		[in]	_numBuildingNormals					- number of buildings's normals
		*	\param		[in]	_numSunNormals					    - number of of suns's normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspotNum(int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals);
		
		/*! \brief This method sets the normals array
		*
		*	\param		[in]	_receiverNormals					- pointer to container of receivers's normals
		*	\param		[in]	_BuildingNormals					- pointer to container of buildings's normals
		*	\param		[in]	_SunNormals					        - pointer to container of suns's normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspotNormals(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals);

		/*! \brief This method sets the positions array
		*
		*	\param		[in]	_receiverPos					    - pointer to container of receivers's positions
		*	\param		[in]	_BuildingPos					    - pointer to container of buildings's positions
		*	\param		[in]	_SunPos					            - pointer to container of suns's positions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspotPosition(const float* _receiverPos, const float* _BuildingPos, const float* _SunPos);

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

		/*! \brief This method sets the new sun vector for hotspot
		*
		*	\param		[in]	newSunVect					     - new sun vector
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspotSun(zVector newSunVect);

		// ------ occlusion -------
		/*! \brief This method sets the new sun vector for hotspot
		*
		*	\param		[in]	_bPositions					     - point to container of buffer positions
		* 	\param		[in]	_bVals					         - number of buffer value
		* 	\param		[in]	_numBuffer					     - number of buffer
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setBuffer(float* _bPositions, int* _bVals, int _numBuffer);

		/*! \brief This method sets the new sun vector for hotspot
		*
		*	\param		[in]	_fPositions					     - point to container of face positions
		* 	\param		[in]	_numFace					     - number of face
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setFacePositions(float* _fPositions, int _numFace);

		//-------------------------------------------------------------------------------
		//---- GET METHODS		
		
		// ------ solar -------
		/*! \brief This method returns the size of container of normals
		*   \returns int  - the size of container of normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numSolarNormals();
		
		/*! \brief This method returns the size of container of sun vectors
		*   \returns int  - the size of container of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numSolarSunVecs();

		/*! \brief This method returns the number of epw data points in the container
		*   \returns int  - the number of epw data points in the container
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numEPWDataPoints();
		
		/*! \brief This method returns the pointer to container of normals
		*   \returns float*  - the pointer to container of normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getSolarRawNormals();
	
		/*! \brief This method returns the pointer to container of colors
		*   \returns float*  - the pointer to container of colors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawColors();
		
		/*! \brief This method returns the pointer to container of struct norm_sunvecs
		*   \returns float*  - the pointer to container of struct norm_sunvecs
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawNormals_SunVectors();
		/*! \brief This method returns the pointer to container of sun vectors for each hour of the year
		*   \returns float*  - the pointer to container of sun vectors for each hour of the year
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawSunVectors_hour();

		/*! \brief This method returns the pointer to container of sun vectors for each hour of 7 special days
		*   \returns float*  - the pointer to container of sun vectors for each hour of 7 special days
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawSunVectors_day();

		/*! \brief This method returns the pointer to container of compass display vectors
		*   \returns float*  - the pointer to container of compass display vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawCompassPts();

		/*! \brief This method returns the pointer to container of EPW data
		*   \returns float*  - the pointer to container of EPW data
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawEPWRadiation();

		/*! \brief This method returns the pointer to container of Cummulative Radiation
		*   \returns float*  - the pointer to container of Cummulative Radiation
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawCummulativeRadiation();

		/*! \brief This method returns the position of sun based on the input date information
		*	\param		[in]	date					    - date information, including year, month, day, hour and minute
		*   \returns zVector  - the position of sun
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zVector getSunPosition(zDate& date);
		
		/*! \brief This method returns the domain of date from sunrise to sunset of one day based on the input date information
		*	\param		[in]	date					    - date information, including year, month, day, hour and minute
		*   \returns zDomainDate  - the domain of date from sunrise to sunset of one day
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zDomainDate getSunRise_SunSet(zDate& date);

		/*! \brief This method returns the domain of dates
		*   \returns zDomainDate  - the domain of dates
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zDomainDate getDomain_Dates();

		/*! \brief This method returns the domain of colors
		*   \returns zDomainColor  - the domain of colors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zDomainColor getDomain_Colors();

		/*! \brief This method returns the geographical location of site, including timeZone, longitude and latitude
		*   \returns zLocation     - the geographical location of site, including timeZone, longitude and latitude
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zLocation getLocation();

		// ------ hotspot -------		
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
		ZSPACE_CUDA_CALLABLE float* getHotspotRawNormals();
				
		/*! \brief This method returns the pointer to container of positions
		*   \returns float*  - pointer to container of positions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getHotspotRawPositions();
						
		/*! \brief This method returns the pointer to container of receivers' scores based on all indexes in the container of sun vectors
		*   \return float*  - pointer to container of receivers' scores based on all indexes in the container of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getAllReceiverScores();

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
		
		/*! \brief This method returns the pointer to container of buildingScore
		*   \returns zBuildingScore*  - the pointer to container of buildingScore
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zBuildingScore* getRawBuildingScore();

		// ------ occlusion -------
		/*! \brief This method returns the number of buffer
		*   \returns int  - the number of buffer
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numBuffer();

		// ------ occlusion -------
		/*! \brief This method returns the number of facePositions
		*   \returns int  - the number of facePositions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numFacePositions();

		/*! \brief This method returns the pointer to container of buffer
		*   \returns float*  - the pointer to container of buffer
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawBuffer();

		/*! \brief This method returns the pointer to container of bufferValues
		*   \returns int*  - the pointer to container of bufferValues
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int* getRawBufferValues();

		/*! \brief This method returns the pointer to container of FacePositions
		*   \returns float*  - the pointer to container of FacePositions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawFacePositions();

		/*! \brief This method returns the pointer to container of Buffer_FacePositions
		*   \returns float*  - the pointer to container of Buffer_FacePositions
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawBuffer_FacePositions();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------	
		
		// ------ solar -------
		/*! \brief This method computes all sun vectors array
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeSunVectors_Year();

		/*! \brief This method computes the sun vectors for each hour of a year
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void computeSunVectors_Hour();

		/*! \brief This method computes the sun vectors for each hour of 7 special days
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void computeSunVectors_Day();

		/*! \brief This method computes the sun vector for an input Date
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST zVector computeSunVector_Date(zDate currentDate);
		
		/*! \brief This method computes the compass display vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeCompass();
		
		/*! \brief This method computes the Cummulative Radiation
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeCummulativeRadiation();

		// ------ hotspot -------	
		/*! \brief This method cleans the previous computed data of Hotspot before another computation
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void cleanScore();

		/*! \brief This method determines whether given three input vectors satisfy the reflection condition within the tolerance range
		*	\param		[in]	incidenceVector					    - input incidence vector, its direction is from sun to building
		*  	\param		[in]	normalVector					    - input normal vector, its direction is same as vertex normal
		*  	\param		[in]	reflectVector					    - input Reflection vector to be check, its direction is from building to receiver
		*  	\param		[in]	angleTolerance					    - angle tolerance with range of [0.0, 1.0], 0.0 means no tolerance
		*   \returns bool  - true if given three vectors satisfy the reflection condition within the tolerance range, else false.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE bool isReflectedVector(zVector incidenceVector, zVector normalVector, zVector reflectVector, float angleTolerance);

		/*! \brief This method computes receivers' scores based on GPU output data
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeReceiverScore();

		
		//--------------------------
		//---- SET MEMORY
		//--------------------------	
	protected:
		
		/*! \brief This method sets memories for containers of solar radiation analysis */
		ZSPACE_CUDA_CALLABLE_HOST void setSolarMemory();

		/*! \brief This method sets memories for containers of Hotspot analysis*/
		ZSPACE_CUDA_CALLABLE_HOST void setHotspotMemory();
				
		/*! \brief This method sets memories for containers of occlusion analysis*/
		ZSPACE_CUDA_CALLABLE_HOST void setOcclusionMemory();

	};
}

#if defined(ZSPACE_CUDA_STATIC_LIBRARY)  || defined(ZSPACE_CUDA_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/energy/zTsSolar.cpp>
#endif

#endif