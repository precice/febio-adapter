//
// Created by fritz on 5/5/22.
//

#include "FEBioReflection.h"
//
// Created by fritz on 5/5/22.
//

#include <FEBioMix/FEBiphasicTPM.h>
#include "FEBioReflection.h"

#include <FEBioMix/FESolutesMaterialPointTPM.h>
#include <FEBioMix/FEMultiphasicMultigeneration.h>
#include <FEBioMix/FEMultiphasicStandardTPM.h>
#include "PreciceCallback.h"
#include <rttr/registration>
#include <FEBioFluid/FEFluidMaterialPoint.h>

RTTR_REGISTRATION {

    registration::class_<FEBiphasicMaterialPointTPM>("FEBiphasicMaterialPointTPM")
            .property("m_phi0", &FEBiphasicMaterialPointTPM::m_phi0);
    registration::class_<FESolutesMaterialPointTPM>("FESolutesMaterialPointTPM")
            .property("m_ca", &FESolutesMaterialPointTPM::m_ca)
            .property("m_psi", &FESolutesMaterialPointTPM::m_psi)
            .method("GetCyp2e1", &FESolutesMaterialPointTPM::GetCyp2e1)
            .method("Vext", &FESolutesMaterialPointTPM::Vext)
            .method("set_concentrations", &FESolutesMaterialPointTPM::set_concentrations)
            .method("set_concentrations_tangents", &FESolutesMaterialPointTPM::set_concentrations_tangents)
            .method("GetViFat", &FESolutesMaterialPointTPM::GetViFat)
            .method("GetViNoFat", &FESolutesMaterialPointTPM::GetViNoFat);
    registration::class_<FEMultiphasicStandardTPM>("FEMultiphasicStandardTPM")
            .method("FluidDensity", &FEMultiphasicStandardTPM::FluidDensity)
            .method("Porosity", &FEMultiphasicStandardTPM::Porosity);
    registration::class_<FEMaterial>("FEMaterial")
            .method("GetID", &FEMaterial::GetID);
    registration::class_<FEFluidMaterialPoint>("FEFluidMaterialPoint")
            .property("m_pf", &FEFluidMaterialPoint::m_pf)
            .property("m_aft", &FEFluidMaterialPoint::m_aft)
            .property("m_vft", &FEFluidMaterialPoint::m_vft);
}
