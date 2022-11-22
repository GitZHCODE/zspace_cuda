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

#ifndef ZSPACE_CD_EXTERN_UTILITIES
#define ZSPACE_CD_EXTERN_UTILITIES

#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cooperative_groups.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_math.h>

#include<headers/zCudaToolsets/base/zCudaDefinitions.h>
#include<headers/zCore/base/zInline.h>
#include<headers/zCore/base/zExtern.h>

#include <stdio.h>
using namespace std;

//---- UTILITIES

ZSPACE_CUDA_EXTERN bool checkCudaExists(string& version)
{

	int runtime;
	cudaRuntimeGetVersion(&runtime);

	cout << "\n CUDA: " << runtime;

	if (runtime >= 10020)
	{
		version = "NVIDIA CUDA 10.2 or higher Installed.";
		return true;
	}
	else
	{
		version = "Install NVIDIA CUDA 10.2 or higher.";
		return false;
	}

	return false;
}



#endif