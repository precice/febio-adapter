// FEBioTPM.cpp : Defines the exported functions for the DLL application.
//

#include <FECore/sdk.h>
#include "PreciceCallback.h"
#include "iostream"

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 1;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
    REGISTER_FECORE_CLASS(PreciceCallback, "precice_callback");
}