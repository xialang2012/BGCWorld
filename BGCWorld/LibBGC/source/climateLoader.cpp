#include "climateLoader.h"

climateLoaderBase::~climateLoaderBase()
{
	if (fileObj != nullptr)
		fileObj->close();
}


void fluxStation::readData()
{

}

fluxStation::~fluxStation()
{
}


StationDataBase::StationDataBase()
{
}

StationDataBase::~StationDataBase()
{
}



double StationDataFlux::calppfd()
{
	double ppfd = 0;
	return ppfd;
}

StationDataFlux::StationDataFlux()
{
}

StationDataFlux::~StationDataFlux()
{
}
