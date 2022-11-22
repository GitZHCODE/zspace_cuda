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

#include <headers/zCudaToolsets/energy/zTsSolar.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSolar::zTsSolar() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSolar::~zTsSolar() {}

	//-------------------------------------------------------------------------------
	//---- SET METHODS

	// ------ solar -------
	ZSPACE_INLINE void zTsSolar::setSolarNormals(const float* _normals, int _numNormals, bool EPWread)
	{
		solarNormals = new float[_numNormals];;
		std::copy(_normals, _normals + _numNormals, solarNormals);

		numSolarNorms = _numNormals;

		setSolarMemory();

		std::copy(solarNormals, solarNormals + (numSolarNorms), norm_sunvecs);
		std::copy(sunVecs_hour, sunVecs_hour + MAX_SUNVECS_HOUR, norm_sunvecs + numSolarNorms);
		//std::copy(sunVecs_date, sunVecs_date + MAX_SUNVECS_DATE, norm_sunvecs + numSolarNorms);

		if (EPWread)
		{
			std::copy(solarNormals, solarNormals + (numSolarNorms), cummulativeRadiation);
			std::copy(epwData_radiation, epwData_radiation + MAX_SUNVECS_HOUR, cummulativeRadiation + numSolarNorms);
		}
	}

	ZSPACE_INLINE bool zTsSolar::setEPWData(string path)
	{
		ifstream myfile;
		myfile.open(path.c_str());

		if (myfile.fail())
		{
			cout << " error in opening file  " << path.c_str() << endl;
			return false;
		}

		epwData_radiation = new float[MAX_SUNVECS_HOUR];

		bool leapYear = (dDate.min.tm_year % 4 == 0) ? true : false;

		bool startCount = false;
		int count = 0;

		while (!myfile.eof())
		{
			string str;
			getline(myfile, str);

			vector<string> perlineData = coreUtils.splitString(str, ",");

			if (perlineData.size() > 0)
			{
				if (perlineData[0] == "LOCATION")
				{
					//location.location = perlineData[1];
					location.latitude = atof(perlineData[6].c_str());
					location.longitude = atof(perlineData[7].c_str());
					location.timeZone = atoi(perlineData[8].c_str());
				}

				if (startCount)
				{

					epwData_radiation[count * 3 + 0] = atof(perlineData[5].c_str()); // temperature
					epwData_radiation[count * 3 + 1] = atof(perlineData[8].c_str()); //pressure
					epwData_radiation[count * 3 + 2] = atof(perlineData[11].c_str()); //radiation

					if (leapYear && count == (59 * 24)) // Feb 28
					{
						for (int k = 0; k < 24; k++)
						{
							count++;

							epwData_radiation[count * 3 + 0] = epwData_radiation[(count - 24) * 3 + 0]; // temperature
							epwData_radiation[count * 3 + 1] = epwData_radiation[(count - 24) * 3 + 1]; //pressure
							epwData_radiation[count * 3 + 2] = epwData_radiation[(count - 24) * 3 + 2]; //radiation
						}

					}

					if (!leapYear && count == (365 * 24)) // Dec 31
					{
						for (int k = 0; k < 24; k++)
						{
							count++;

							epwData_radiation[count * 3 + 0] = INVALID_VAL; // temperature
							epwData_radiation[count * 3 + 1] = INVALID_VAL; //pressure
							epwData_radiation[count * 3 + 2] = INVALID_VAL; //radiation
						}
					}

					//printf("\n %i | %1.2f ", count * 3 + 2, epwData_radiation[count * 3 + 2]);

					//epwData[count].humidity = atof(perlineData[7].c_str());
					//epwData[count].windDirection = atof(perlineData[19].c_str());
					//epwData[count].windSpeed = atof(perlineData[20].c_str());							

					count++;
				}

				if (perlineData[0] == "DATA PERIODS") startCount = true;
			}

		}


		printf("\n count : %i ", count);
		myfile.close();

		return true;
	}

	ZSPACE_INLINE void zTsSolar::setDomain_Dates(zDomainDate& _dDate)
	{
		dDate = _dDate;
	}

	ZSPACE_INLINE void zTsSolar::setDomain_Colors(zDomainColor& _dColor)
	{
		dColor = _dColor;
	}

	ZSPACE_INLINE void zTsSolar::setLocation(zLocation& _location)
	{
		location = _location;
	}

	// ------ hotspot -------
	ZSPACE_INLINE void zTsSolar::setHotspotNum(int _numReceiverNormals, int _numBuildingNormals, int _numSunNormals)
	{
		numBuildingNormals = _numBuildingNormals;
		numSunNormals = _numSunNormals;
		numReceiverNormals = _numReceiverNormals;
		setHotspotMemory();
	}

	ZSPACE_INLINE void zTsSolar::setHotspotNormals
	(const float* _receiverNormals, const float* _BuildingNormals, const float* _SunNormals)
	{
		std::copy(_BuildingNormals, _BuildingNormals + numBuildingNormals, hotspotNormals);
		std::copy(_SunNormals, _SunNormals + numSunNormals, hotspotNormals + numBuildingNormals);
		std::copy(_receiverNormals, _receiverNormals + numReceiverNormals, hotspotNormals + numSunNormals + numBuildingNormals);
	}

	ZSPACE_INLINE void zTsSolar::setHotspotPosition
	(const float* _receiverPos, const float* _BuildingPos, const float* _SunPos)
	{
		std::copy(_BuildingPos, _BuildingPos + numBuildingNormals, hotspotPositions);
		std::copy(_SunPos, _SunPos + numSunNormals, hotspotPositions + numBuildingNormals);
		std::copy(_receiverPos, _receiverPos + numReceiverNormals, hotspotPositions + numSunNormals + numBuildingNormals);
	}

	ZSPACE_INLINE void zTsSolar::setReceiverScores(const float* _receiverScores)
	{
		std::copy(_receiverScores, _receiverScores + (numReceiverNormals / 3), receiverScores);
	}

	ZSPACE_INLINE void zTsSolar::setOcclusion(const float* _occlusion)
	{
		std::copy(_occlusion, _occlusion + (numBuildingNormals / 3), occlusion);
	}

	ZSPACE_INLINE void zTsSolar::setHotspotSun(zVector newSunVect)
	{
		hotspotNormals[numBuildingNormals] = newSunVect.x;
		hotspotNormals[numBuildingNormals+1] = newSunVect.y;
		hotspotNormals[numBuildingNormals+2] = newSunVect.z;
	}

	// ------ occlusion -------
	ZSPACE_INLINE void zTsSolar::setBuffer(float* _bPositions, int* _bVals, int _numBuffer)
	{
		bufferPositions = new float[_numBuffer];;
		std::copy(_bPositions, _bPositions + (_numBuffer), bufferPositions);

		bufferValues = new int[_numBuffer];
		std::copy(_bVals, _bVals + (_numBuffer), bufferValues);

		num_Buffer = _numBuffer;
	}

	ZSPACE_INLINE void zTsSolar::setFacePositions(float* _fPositions, int _numFace)
	{
		facePositions = new float[_numFace];;
		std::copy(_fPositions, _fPositions + _numFace, facePositions);

		numFaces = _numFace;

		setOcclusionMemory();

		std::copy(bufferPositions, bufferPositions + (num_Buffer), buffer_facePositions);
		std::copy(facePositions, facePositions + numFaces, buffer_facePositions + num_Buffer);
	}

	//-------------------------------------------------------------------------------
	//---- GET METHODS	
	
	// ------ solar -------
	ZSPACE_INLINE int zTsSolar::numSolarNormals()
	{
		return numSolarNorms;
	}

	ZSPACE_INLINE int zTsSolar::numSolarSunVecs()
	{
		return MAX_SUNVECS_HOUR;
	}

	ZSPACE_INLINE int zTsSolar::numEPWDataPoints()
	{
		return numEPWData;
	}

	ZSPACE_INLINE float* zTsSolar::getSolarRawNormals()
	{
		return solarNormals;
	}

	ZSPACE_INLINE float* zTsSolar::getRawColors()
	{
		return colors;
	}

	ZSPACE_INLINE float* zTsSolar::getRawNormals_SunVectors()
	{
		return norm_sunvecs;
	}

	ZSPACE_INLINE float* zTsSolar::getRawSunVectors_hour()
	{
		return sunVecs_hour;
	}

	ZSPACE_INLINE float* zTsSolar::getRawSunVectors_day()
	{
		return sunVecs_days;
	}

	ZSPACE_INLINE float* zTsSolar::getRawCompassPts()
	{
		return compassPts;
	}

	ZSPACE_INLINE float* zTsSolar::getRawEPWRadiation()
	{
		return epwData_radiation;
	}

	ZSPACE_INLINE float* zTsSolar::getRawCummulativeRadiation()
	{
		return cummulativeRadiation;
	}

	ZSPACE_INLINE zVector zTsSolar::getSunPosition(zDate& date)
	{
		float LocalTime = date.tm_hour + (date.tm_min / 60.0);

		double JD = date.toJulian();

		double phi = location.latitude;
		double lambda = location.longitude;

		double n = JD - 2451545.0;

		double LDeg = (double)fmod((280.460 + 0.9856474 * n), 360.0);
		double gDeg = (double)fmod((357.528 + 0.9856003 * n), 360.0);

		double LambdaDeg = LDeg + 1.915 * sin(gDeg * DEG_TO_RAD) + 0.01997 * sin(2 * gDeg * DEG_TO_RAD);

		double epsilonDeg = 23.439 - 0.0000004 * n;

		double alphaDeg;
		alphaDeg = atan(cos(epsilonDeg * DEG_TO_RAD) * tan(LambdaDeg * DEG_TO_RAD));
		alphaDeg *= RAD_TO_DEG;
		if (cos(LambdaDeg * DEG_TO_RAD) < 0)	alphaDeg += (4 * (atan(1.0) * RAD_TO_DEG));

		double deltaDeg = asin(sin(epsilonDeg * DEG_TO_RAD) * sin(LambdaDeg * DEG_TO_RAD)) * RAD_TO_DEG;

		zDate dZero(date.tm_year, date.tm_mon, date.tm_mday, 0, 0);
		double JDNull = dZero.toJulian();

		double TNull = ((JDNull - 2451545.0) / 36525);
		double T = LocalTime - location.timeZone;

		double thetaGh = 6.697376 + 2400.05134 * TNull + 1.002738 * T;

		double thetaG = (double)fmod(thetaGh * 15.0, 360.0);
		double theta = thetaG + lambda;

		double tauDeg = theta - alphaDeg;

		double denom = (cos(tauDeg * DEG_TO_RAD) * sin(phi * DEG_TO_RAD) - tan(deltaDeg * DEG_TO_RAD) * cos(phi * DEG_TO_RAD));
		double aDeg = atan(sin(tauDeg * DEG_TO_RAD) / denom);
		aDeg *= RAD_TO_DEG;
		if (denom < 0) aDeg = aDeg + 180;
		aDeg += 180; //add 180 to azimuth to compute from the north.

		double hDeg = asin(cos(deltaDeg * DEG_TO_RAD) * cos(tauDeg * DEG_TO_RAD) * cos(phi * DEG_TO_RAD) + sin(deltaDeg * DEG_TO_RAD) * sin(phi * DEG_TO_RAD));
		hDeg *= RAD_TO_DEG;

		double valDeg = hDeg + (10.3 / (hDeg + 5.11));
		double RDeg = 1.02 / (tan(valDeg * DEG_TO_RAD));

		double hRDeg = hDeg + (RDeg / 60);

		return coreUtils.sphericalToCartesian(aDeg, hRDeg, 1.0);
	}

	ZSPACE_INLINE zDomainDate zTsSolar::getSunRise_SunSet(zDate& date)
	{
		zDomainDate out;

		zDate temp = date;
		temp.tm_hour = 12;
		temp.tm_min = 0;
		temp.tm_sec = 0;

		double jd = temp.toJulian();

		double n = jd - 2451545.0 + 0.0008;

		double js = n - (location.longitude / 360.0);

		double m = (double)fmod((357.5291 + 0.98560028 * js), 360.0) * DEG_TO_RAD; // radians

		double c = 1.9148 * sin(m) + 0.02 * sin(2 * m) + 0.0003 * sin(3 * m);

		double mDeg = m * RAD_TO_DEG;
		double lambdaDeg = (double)fmod((mDeg + c + 180 + 102.9372), 360.0); //deg
		double lambda = lambdaDeg * DEG_TO_RAD;

		double jt = 2451545.0 + js + ((0.0053 * sin(m)) - (0.0069 * sin(2 * lambda)));

		double delta = asin(sin(lambda) * sin(23.44 * DEG_TO_RAD));

		double cosOmega = (sin(-0.83 * DEG_TO_RAD) - (sin(location.latitude * DEG_TO_RAD) * sin(delta))) / (cos(location.latitude * DEG_TO_RAD) * cos(delta));
		double omegaDeg = acos(cosOmega) * RAD_TO_DEG;

		double j;

		//sunrise
		j = jt - (omegaDeg / 360.0);
		out.min.fromJulian(j);

		//sunset
		j = jt + (omegaDeg / 360.0);
		out.max.fromJulian(j);
		return out;
	}

	ZSPACE_INLINE zDomainDate zTsSolar::getDomain_Dates()
	{
		return dDate;
	}

	ZSPACE_INLINE zDomainColor zTsSolar::getDomain_Colors()
	{
		return dColor;
	}

	ZSPACE_INLINE zLocation zTsSolar::getLocation()
	{
		return location;
	}

	// ------ hotspot -------
	ZSPACE_INLINE int zTsSolar::getNumSunNormals()
	{
		return numSunNormals;
	}

	ZSPACE_INLINE int zTsSolar::getNumBuildingNormals()
	{
		return numBuildingNormals;
	}

	ZSPACE_INLINE int zTsSolar::getNumReceiverNormals()
	{
		return numReceiverNormals;
	}

	ZSPACE_INLINE float* zTsSolar::getHotspotRawNormals()
	{
		return hotspotNormals;
	}

	ZSPACE_INLINE float* zTsSolar::getHotspotRawPositions()
	{
		return hotspotPositions;
	}

	ZSPACE_INLINE float* zTsSolar::getAllReceiverScores()
	{
		return receiverScores;
	}

	ZSPACE_INLINE float zTsSolar::getAllHighScore()
	{
		float highestScore = 0;

		for (int i = 0; i < (numBuildingNormals / 3); i += 1)
		{
			if (receiverScores[i] > highestScore)
			{
				highestScore = receiverScores[i];
			}
		}

		return highestScore;
	}

	ZSPACE_INLINE float zTsSolar::getAllScoreSum()
	{
		float sum = 0;

		for (int i = 0; i < (numBuildingNormals / 3); i += 1)
		{
			sum += receiverScores[numReceiverNormals + i];
		}
		
		return sum;
	}

	ZSPACE_INLINE float* zTsSolar::getRawOcclusion()
	{
		return occlusion;
	}

	ZSPACE_INLINE zBuildingScore* zTsSolar::getRawBuildingScore()
	{
		return buildingScore;
	}
	
	// ------ occlusion -------
	ZSPACE_INLINE int zTsSolar::numBuffer()
	{
		return num_Buffer;
	}

	ZSPACE_INLINE int zTsSolar::numFacePositions()
	{
		return numFaces;
	}

	ZSPACE_INLINE float* zTsSolar::getRawBuffer()
	{
		return bufferPositions;
	}

	ZSPACE_INLINE int* zTsSolar::getRawBufferValues()
	{
		return bufferValues;
	}

	ZSPACE_INLINE float* zTsSolar::getRawFacePositions()
	{
		return facePositions;
	}

	ZSPACE_INLINE float* zTsSolar::getRawBuffer_FacePositions()
	{
		return buffer_facePositions;
	}
	
	//---- COMPUTE METHODS	

	// ------ solar -------
    ZSPACE_INLINE void zTsSolar::computeSunVectors_Year()
	{
		computeSunVectors_Hour();
		computeSunVectors_Day();
	}

	ZSPACE_INLINE void zTsSolar::computeSunVectors_Hour()
	{
		zDate min = dDate.min;
		min.tm_mon = 1;
		min.tm_mday = 1;
		min.tm_hour = 0;
		min.tm_min = 0;


		zDate max = dDate.max;
		max.tm_mon = 12;
		max.tm_mday = 31;
		max.tm_hour = 23;
		max.tm_min = 59;


		time_t  unixTime_s = min.toUnix();
		time_t  unixTime_e = max.toUnix();

		// get minute domain per day
		zDate minHour = min;
		zDate maxHour(min.tm_year, min.tm_mon, min.tm_mday, max.tm_hour, max.tm_min);

		time_t  unixTime_sh = minHour.toUnix();
		time_t  unixTime_eh = maxHour.toUnix();


		//get total number of vectors
		sunVecs_hour = new float[MAX_SUNVECS_HOUR];

		int hrCount = 0;

		for (time_t hour = unixTime_sh; hour <= unixTime_eh; hour += 3600)
		{
			int dCount = 0;
			for (time_t day = unixTime_s; day <= unixTime_e; day += 86400)
			{
				zDate currentDate;
				currentDate.fromUnix(day + hour - unixTime_s);;

				zVector sunPos = getSunPosition(currentDate);

				if (sunPos.z >= 0)
				{
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 0] = sunPos.x;
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 1] = sunPos.y;
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 2] = sunPos.z;
				}
				else
				{
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 0] = INVALID_VAL;
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 1] = INVALID_VAL;
					sunVecs_hour[(((hrCount * 366) + dCount) * 3) + 2] = INVALID_VAL;
				}

				dCount++; ;
			}

			hrCount++;
		}

		// for non leap years
		if (min.tm_year % 4 != 0)
		{
			for (int i = (365 * 24 * 3); i < MAX_SUNVECS_HOUR; i++)
			{
				sunVecs_hour[i] = INVALID_VAL;
			}
		}
	}

	ZSPACE_INLINE void zTsSolar::computeSunVectors_Day()
	{

		zDate days[7];

		int _year = dDate.min.tm_year;

		days[0] = zDate(_year, 06, 21, 0, 1);
		days[1] = zDate(_year, 12, 21, 0, 1);
		days[2] = zDate(_year, 1, 28, 0, 1);
		days[3] = zDate(_year, 2, 28, 0, 1);
		days[4] = zDate(_year, 3, 21, 0, 1);
		days[5] = zDate(_year, 4, 15, 0, 1);
		days[6] = zDate(_year, 5, 15, 0, 1);

		sunVecs_days = new float[MAX_SUNVECS_DAY];

		for (int i = 0; i < 7; i++)
		{

			zDate init = days[i];
			zDate end = days[i];

			init.tm_hour = 0;
			init.tm_min = 1;
			end.tm_hour = 23;
			end.tm_min = 59;

			time_t  unixTime_s = init.toUnix();
			time_t  unixTime_e = end.toUnix();

			int count = 0;
			for (time_t hour = unixTime_s; hour <= unixTime_e; hour += 3600)
			{
				zDate currentDate;
				currentDate.fromUnix(hour);
				zDomainDate dd = getSunRise_SunSet(currentDate);


				zVector sunPos;

				sunPos = getSunPosition(currentDate);
				if (sunPos.z >= 0)
				{
					sunVecs_days[(((i * 24) + count) * 3) + 0] = sunPos.x;
					sunVecs_days[(((i * 24) + count) * 3) + 1] = sunPos.y;
					sunVecs_days[(((i * 24) + count) * 3) + 2] = sunPos.z;
				}
				else
				{
					sunVecs_days[(((i * 24) + count) * 3) + 0] = INVALID_VAL;
					sunVecs_days[(((i * 24) + count) * 3) + 1] = INVALID_VAL;
					sunVecs_days[(((i * 24) + count) * 3) + 2] = INVALID_VAL;
				}


				if (currentDate.tm_hour == dd.min.tm_hour)
				{
					sunPos = getSunPosition(dd.min);
					sunVecs_days[(((i * 24) + count) * 3) + 0] = sunPos.x;
					sunVecs_days[(((i * 24) + count) * 3) + 1] = sunPos.y;
					sunVecs_days[(((i * 24) + count) * 3) + 2] = sunPos.z;
				}

				if (currentDate.tm_hour == dd.max.tm_hour + 1)
				{
					sunPos = getSunPosition(dd.max);
					sunVecs_days[(((i * 24) + count) * 3) + 0] = sunPos.x;
					sunVecs_days[(((i * 24) + count) * 3) + 1] = sunPos.y;
					sunVecs_days[(((i * 24) + count) * 3) + 2] = sunPos.z;
				}

				count++;
			}
		}
	}

	ZSPACE_INLINE zVector zTsSolar::computeSunVector_Date(zDate currentDate)
	{
		sunVecs_date = new float[MAX_SUNVECS_DATE];
		zVector sunPos;
		sunPos = getSunPosition(currentDate);

		if (sunPos.z < 0)
		{
			sunPos.x = INVALID_VAL;
			sunPos.y = INVALID_VAL;
			sunPos.z = INVALID_VAL;
		}

		sunVecs_date[0] = sunPos.x;
		sunVecs_date[1] = sunPos.y;
		sunVecs_date[2] = sunPos.z;
		
		sunPos *= -1;
		return sunPos;
	}

	ZSPACE_INLINE void zTsSolar::computeCompass()
	{
		compassPts = new float[COMPASS_SUBD];
		float deg = (float)(360 / (float)(12));

		//zVector pos(0, 1, 0);

		for (int i = 0; i < 2; i++)
		{
			//if (i > 0) pos *= 1.1;
			zVector pos(0, 1 + (i * 0.1), 0);

			for (int j = 0; j < 12; j++)
			{
				int id = (i * 12 + j) * 3;
				compassPts[id] = pos.x;
				compassPts[id + 1] = pos.y;
				compassPts[id + 2] = pos.z;
				pos = pos.rotateAboutAxis(zVector(0, 0, 1), deg);
			}

		}
	}

	ZSPACE_INLINE void zTsSolar::computeCummulativeRadiation()
	{
		for (int i = 0; i < numSolarNorms; i += 3)

		{
			zVector norm(norm_sunvecs[i + 0], norm_sunvecs[i + 1], norm_sunvecs[i + 2]);

			int sunvecs_offset = numSolarNorms - i;

			float angle = 0;
			int count = 0;
			for (int o = i; o < i + MAX_SUNVECS_HOUR; o += 3)
			{
				int j = o + sunvecs_offset;
				zVector sVec(norm_sunvecs[j + 0], norm_sunvecs[j + 1], norm_sunvecs[j + 2]);

				if (sVec.x != INVALID_VAL && sVec.y != INVALID_VAL && sVec.z != INVALID_VAL)
				{
					angle += norm.angle(sVec);
					count++;
				}


			}
			angle /= count;


			if (angle > 90.0)
			{
				colors[i + 0] = dColor.min.h;
				colors[i + 1] = dColor.min.s;
				colors[i + 2] = dColor.min.v;
			}
			else
			{
				colors[i + 0] = coreUtils.ofMap(angle, 90.0f, 0.0f, dColor.min.h, dColor.max.h);
				colors[i + 1] = coreUtils.ofMap(angle, 90.0f, 0.0f, dColor.min.s, dColor.max.s);
				colors[i + 2] = coreUtils.ofMap(angle, 90.0f, 0.0f, dColor.min.v, dColor.max.v);
			}



		}
	}

	// ------ hotspot -------	
	ZSPACE_INLINE void zTsSolar::cleanScore()
	{
		for (int j = 0; j < (numReceiverNormals/3); j += 1)
		{
			receiverScores[j] = 0;
		}

		for (int i = 0; i < numBuildingNormals; i += 3)
		{
			int buildID = floor(i / 3);
			buildingScore[buildID].receiverCount = 0;
		}

	}

	ZSPACE_INLINE bool zTsSolar::isReflectedVector(zVector incidenceVector, zVector normalVector, zVector reflectVector, float angleTolerance)
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

	ZSPACE_INLINE void zTsSolar::computeReceiverScore()
	{
		int maxCount = 0;

		for (int i = 0; i < numBuildingNormals; i += 3)
		{
			int buildID = floor(i / 3);

			for (int r = 0; r < buildingScore[buildID].receiverCount; r++)
			{
				int receiverID = buildingScore[buildID].receiverID[r];
				receiverScores[receiverID]++;
			}

			if (maxCount < buildingScore[buildID].receiverCount)
			{
				maxCount = buildingScore[buildID].receiverCount;
			}
		}

		printf("\n maxCount of receiver for each buildingVertex = %i", maxCount);
	}

	//---- PROTECTED METHODS
	ZSPACE_INLINE void zTsSolar::setSolarMemory()
	{
		if ((numSolarNorms + MAX_SUNVECS_HOUR) < solarMemSize) return;
		else
		{
			while (solarMemSize < (numSolarNorms + MAX_SUNVECS_HOUR)) solarMemSize += d_MEMORYMULTIPLIER;

			norm_sunvecs = new float[solarMemSize];
			cummulativeRadiation = new float[solarMemSize];

			// set to  Num Normals
			colors = new float[solarMemSize];
		}
	}

	ZSPACE_INLINE void zTsSolar::setHotspotMemory()
	{
		if ((numBuildingNormals + numSunNormals + numReceiverNormals) < hotspotMemSize) return;
		else
		{
			while (hotspotMemSize < numBuildingNormals + numSunNormals + numReceiverNormals) hotspotMemSize += d_MEMORYMULTIPLIER;
		}

		hotspotNormals = new float[hotspotMemSize];
		hotspotPositions = new float[hotspotMemSize]; 
		occlusion = new float[numBuildingNormals / 3];
		receiverScores = new float[numReceiverNormals / 3];
		buildingScore = new zBuildingScore[numBuildingNormals / 3];
	}

	ZSPACE_INLINE void zTsSolar::setOcclusionMemory()
	{
		if ((num_Buffer + numFaces) < occlusionMemSize) return;
		else
		{
			while (occlusionMemSize < (num_Buffer + numFaces)) occlusionMemSize += d_MEMORYMULTIPLIER;

			buffer_facePositions = new float[occlusionMemSize];
		}
	}
}