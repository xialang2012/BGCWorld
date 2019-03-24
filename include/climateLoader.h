#ifndef CLIMATELOADER_H
#define CLIMATELOADER_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

class climateLoaderBase
{
public:
	climateLoaderBase(const std::string & s) : inFile(s), fileObj(nullptr) {};
	~climateLoaderBase();
	//virtual void OpenData();

	std::ifstream* fileObj;
private:
	std::string inFile;
};

class fluxStation : public climateLoaderBase
{
public:
	fluxStation(const std::string & s) : climateLoaderBase(s) {};
	void readData();

	~fluxStation();

private:

};

class StationDataBase
{
public:
	StationDataBase();
	~StationDataBase();

	int year, month, day, hour, minute, days;   // time flag
	int stationId;  // station id	
	double tempAvg, tempHigh, tempLow;  // temp
	double precipitation; // precipitation
	double windSpeed;  // wind speed

private:

};

// flux station data (half hour data)
class StationDataFlux : StationDataBase
{
public:
	StationDataFlux();
	~StationDataFlux();

	int YEAR, DAY, HRMIN;  //time year, day of year, minutes
	double Rd;  //(umol / m2 / s)
	double GPP;	//(umol / m2 / s)
	double NEE;	//(umol / m2 / s)
	double TA;	//(deg C)
	double PREC;  //(mm)
	double RH;	//(%)
	double VPD;	//(kPa)
	double Srad;  // (W/m2)
	double SWC;	//(%)
	double PAR;	//(umol / m2 / s)

private:
	double calppfd();
};

#endif // !CLIMATELOADER