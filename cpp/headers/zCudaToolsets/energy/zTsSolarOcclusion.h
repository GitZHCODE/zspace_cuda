// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Federico Borello <federico.borello@zaha-hadid.com>, Cesar Fragachan <cesar.fragachan@zaha-hadid.com>
//

#ifndef ZSPACE_TS_SOLAR_OCCLUSION
#define ZSPACE_TS_SOLAR_OCCLUSION

#pragma once
#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zDomain.h>
#include<headers/zCore/utilities/zUtilsCore.h>
#include <headers/zApp/include/zObjects.h>
#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnMeshField.h>
#include <headers/zCore/base/zExtern.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

namespace zSpace
{
	//struct zNorm_SunVec
	//{
	//	zVector norm;
	//	zVector sunVec;
	//};

	class ZSPACE_CUDA zTsSolarOcclusion

	{
	protected:
		//sun variables
		zVector* sunVector;
		zVector* sunCenter;
		zTransform transform;
		zTransform transformLocal;

		//building variables
		zObjMesh buildingMesh;
		vector<size_t> faceID;
		int numVertices;
		int numFaces;

		//field variables
		zObjMeshScalarField oMesh_Field;
		zPoint fieldMinBB, fieldMaxBB;

		double unitX;
		double unitY;
		int resolutionX;
		int resolutionY;

		//output
		zBoolArray vertexOcclusion;
		zIntArray occlusionVal;


		//	//--------------------------
		//	//---- COMPUTE METHODS
		//	//--------------------------
		
		//global
		void init();

		//sun
		void computeTransform();

		//building
		void flatternBuilding();
		void computeFieldBBox();
		void sortFaceByDist();
		void computeFaceProjection(int currentfaceID);
		void computeVertexValue(int currentfaceID);

		//field
		void createField();
		void updateField(int currentfaceID);


	public:

		zUtilsCore core;

		//	/*!	\brief core utilities Object  */
			//--------------------------
			//---- CONSTRUCTOR
			//--------------------------

			/*! \brief Default constructor.
			*	\since version 0.0.4
			*/
		zTsSolarOcclusion();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*	\since version 0.0.4
		*/
		~zTsSolarOcclusion();

		//	//--------------------------
		//	//---- SET METHODS
		//	//--------------------------	
		//		   
		void setSunVec(zVector& _sunVec);
		void setSunCenter(zVector& _sunCenter);
		void setBuildingMesh(zObjMesh _buildingMesh);
		void setResolution(int& _resolutionX, int& _resolutionY);

		//	//--------------------------
		//	//---- GET METHODS
		//	//--------------------------	

		zBoolArray getOcclusion();
		zIntArray getOcclusionValue();
		zPointArray getBBox();
		zObjMeshScalarField getField();
		zPointArray getFace();
		zTransform getTransform();
		zTransform getTransformLocal();



		//	//--------------------------
		//	//---- COMPUTE METHODS
		//	//--------------------------
		// 
		//void computeBBox();

		void computeOcclusion();
	};
}

#if defined(ZSPACE_CUDA_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zCudaToolsets/energy/zTsSolarOcclusion.cpp>
#endif

#endif