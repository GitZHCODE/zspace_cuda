// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2022 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Federico Borello <federico.borello@zaha-hadid.com>, Cesar Fragachan <cesar.fragachan@zaha-hadid.com>, Ling Mao <Ling.Mao@zaha-hadid.com>
//

#ifndef ZSPACE_TS_SOLAR_ANALYSIS
#define ZSPACE_TS_SOLAR_ANALYSIS

#pragma once

#define MAX_SUNVECS_HOUR 366*24 * 3
#define MAX_SUNVECS_DAY 24*7 * 3
#define INVALID_VAL -10000.0

#define COMPASS_SUBD 12*2*3

#define RAD_ANGLE_MIN 30.0
#define RAD_ANGLE_MAX 120.0

#define RAD_MIN 300.0
#define RAD_MAX 1000.0


#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zDomain.h>
#include<headers/zCore/utilities/zUtilsCore.h>

namespace zSpace
{

	struct zNorm_SunVec
	{
		zVector norm;
		zVector sunVec;
	};

	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsSolarAnalysis
	*	\brief A tool set to do solar analysis.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_CUDA_TOOLSET zTsSolarAnalysis
	{
	private:


		/*!	\brief pointer to container of normals	*/
		float *normals;

		/*!	\brief number of face normals in the container	*/
		int numNorms;

		/*!	\brief pointer to container of sun vectors	*/
		float *sunVecs_hour;

		/*!	\brief pointer to container of sun vectors	*/
		float *sunVecs_days;

		/*!	\brief pointer to container of compass display vectors	*/
		float *compassPts;

		/*!	\brief pointer to container of epw data	*/
		float* epwData_radiation;		

		/*!	\brief number of epw data points in the container	*/
		int numData;

		/*!	\brief date domain	*/
		zDomainDate dDate;

		/*!	\brief location information */
		zLocation location;

		/*!	\brief color domain	*/
		zDomainColor dColor;

		/*!	\brief pointer to container of cummulative raditation*/
		float *cummulativeRadiation;

		/*!	\brief pointer to container of colors*/
		float *colors;

		/*!	\brief pointer to container of normals + sunvectors*/
		float *norm_sunvecs;

		/*!	\brief size of container for normals, cummulative raditation, colors*/
		int memSize;


	public:

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;		


		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST zTsSolarAnalysis();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST ~zTsSolarAnalysis();

		//--------------------------
		//---- SET METHODS
		//--------------------------	
			   
		/*! \brief This method sets container of normals
		*	\param		[in]	_normals					- pointer to container of normals
		*	\param		[in]	_numNormals					- number of face normals in the container
		*	\param		[in]	EPWread					    - true if reads EPW as input, else false  
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setNormals(const float *_normals, int _numNormals, bool EPWread);

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
		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Dates(zDomainDate & _dDate);

		/*! \brief This method sets the domain of colors
		*	\param		[in]	_dColor					    - domain of colors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Colors(zDomainColor & _dColor);

		/*! \brief This method sets the geographical location of site, including timeZone, longitude and latitude
		*	\param		[in]	_location					- geographical location, including timeZone, longitude and latitude
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setLocation(zLocation &_location);

		//--------------------------
		//---- GET METHODS
		//--------------------------	

		/*! \brief This method returns the size of container of normals
		*   \returns int  - the size of container of normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numNormals();

		/*! \brief This method returns the size of container of sun vectors
		*   \returns int  - the size of container of sun vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numSunVecs();

		/*! \brief This method returns the number of epw data points in the container
		*   \returns int  - the number of epw data points in the container
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE int numDataPoints();

		/*! \brief This method returns the pointer to container of normals
		*   \returns float*  - the pointer to container of normals
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE float* getRawNormals();

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
		ZSPACE_CUDA_CALLABLE zVector getSunPosition(zDate &date);

		/*! \brief This method returns the domain of date from sunrise to sunset of one day based on the input date information
		*	\param		[in]	date					    - date information, including year, month, day, hour and minute
		*   \returns zDomainDate  - the domain of date from sunrise to sunset of one day
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zDomainDate getSunRise_SunSet(zDate &date);

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

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------	

		/*! \brief This method computes sun vectors of the whole year
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeSunVectors_Year();

		/*! \brief This method computes the compass display vectors
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeCompass();

		/*! \brief This method computes the Cummulative Radiation
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE void computeCummulativeRadiation();

		/*! \brief This method computes and returns the sun vector based on the input date information
		*	\param		[in]	currentDate					- current date information
		*   \returns zVector  - the sun vector based on the input date information
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE zVector computeSunVectors_Hotspot(zDate currentDate);

		//--------------------------
		//---- DISPLAY METHODS
		//--------------------------	

		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------	
	protected:

		/*! \brief This method sets memories for containers in zTsSolarAnalysis
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void setMemory();

		/*! \brief This method computes the sun vectors for each hour of the year
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void computeSunVectors_Hour();

		/*! \brief This method computes the sun vectors for each hour of 7 special days
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST void computeSunVectors_Day();		

	};
}

#if defined(ZSPACE_CUDA_STATIC_LIBRARY)  || defined(ZSPACE_CUDA_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/energy/zTsSolarAnalysis.cpp>
#endif

#endif