// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>,
//

#ifndef ZSPACE_TS_MESHFIELD
#define ZSPACE_TS_MESHFIELD

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

	
	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsMeshField
	*	\brief A tool set for mesh fields.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_CUDA zTsPixelField
	{
	private:

		/*!	\brief pointer to container of vertex positions	*/
		float* vPositions;

		/*!	\brief pointer to container of vertex colors*/
		float* vColors;

		/*!	\brief number of points in the container	*/
		int numV;

		/*!	\brief pointer to container of polygon vertex connects. Quads only	*/
		int* polyVConnects;
		
		/*!	\brief number of points in the container	*/
		int numP;

		/*!	\brief pointer to container of isovalues. */
		float* pixelValues;

		/*!	\brief number of points in the container	*/
		int numE;

		/*!	\brief pointer to container of polygon edge connects. Quads only	*/
		int* polyEConnects;

		/*!	\brief pointer to container of edge colors*/
		float* eColors;

		/*!	\brief pointer to container of contour positions, every pair of 2 positions make an graph edge*/
		float* contourPositions;

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
		ZSPACE_CUDA_CALLABLE_HOST zTsPixelField();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		ZSPACE_CUDA_CALLABLE_HOST ~zTsPixelField();

		//--------------------------
		//---- SET METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE_HOST void setPixels(const float* _vPositions, const float* _vColors, float* _pixelValues, const int* _polyconnects, int _numVertices, int _numPolygons);

		ZSPACE_CUDA_CALLABLE_HOST void setPixelEdges(const int* edgeconnects, const float* _edgeColors, int _numEdges);
		
		ZSPACE_CUDA_CALLABLE_HOST void setDomain_Colors(zDomainColor& _dColor);
				
		//--------------------------
		//---- GET METHODS
		//--------------------------	

		ZSPACE_CUDA_CALLABLE int numVertices();

		ZSPACE_CUDA_CALLABLE int numPolygons();

		ZSPACE_CUDA_CALLABLE int numEdges();

		ZSPACE_CUDA_CALLABLE float* getRawVertexPositions();

		ZSPACE_CUDA_CALLABLE float* getRawVertexColors();

		ZSPACE_CUDA_CALLABLE float* getRawPolyConnects();

		ZSPACE_CUDA_CALLABLE float* getRawEdgeConnects();

		ZSPACE_CUDA_CALLABLE float* getRawEdgeColors();

		ZSPACE_CUDA_CALLABLE float* getRawContours();

		ZSPACE_CUDA_CALLABLE zDomainColor getDomain_Colors();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------	

		//--------------------------
		//---- DISPLAY METHODS
		//--------------------------	

		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------	
	protected:

		ZSPACE_CUDA_CALLABLE_HOST void setMemory();
				

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/field/zTsPixelField.cpp>
#endif

#endif