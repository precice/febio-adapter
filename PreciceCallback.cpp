//
// Created by fritz on 2/25/22.
//

#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FEBioMix/FEBiphasicTPM.h>
#include <FEBioMix/FESolutesMaterialPointTPM.h>
#include <FEBioMix/FEMultiphasicMultigeneration.h>
#include <FECore/log.h>
#include "PreciceCallback.h"
#include "iostream"
#include <nlohmann/json.hpp>
#include <fstream>
#include <variant>
#include <rttr/registration>
#include <FEBioFluid/FEFluidMaterialPoint.h>

using namespace rttr;


std::pair<int, vector<double>>
PreciceCallback::getRelevantMaterialPoints(FEModel *fem, const std::string &elementName) {
    vector<double> vertexPositions;
    int counter = 0;
    FEElementSet *elementSet = fem->GetMesh().FindElementSet(elementName);
    if (!elementSet) {
        feLogError((elementName + std::string("ElementSet not found")).c_str());
        throw FEException("ElementSet not found");
    }
    for (int i = 0; i < elementSet->Elements(); i++) {
        FEElement &element = elementSet->Element(i);
        for (int j = 0; j < element.GaussPoints(); j++) {
            FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
            vec3d coord = materialPoint->m_r0; // We are using the inital position to identify each node, as we are not able to support dynamic meshes yet
            vertexPositions.push_back(coord.x);
            vertexPositions.push_back(coord.y);
            vertexPositions.push_back(coord.z);
            counter++;
        }
    }
    return std::pair(counter, vertexPositions);

}


void PreciceCallback::Init(FEModel *fem) {
    feLogInfo("Init precice");
    using json = nlohmann::json;
    //TODO make this config path configurable
    std::ifstream config_file("febio-config.json");
    json config_json;
    config_file >> config_json;
    precice = new precice::SolverInterface(config_json["coupling_params"]["participant_name"],
                                           config_json["coupling_params"]["config_file_name"], 0, 1);
    dimensions = precice->getDimensions();
    FEMesh &femMesh = fem->GetMesh();
    // Counting all Material Points in the Mesh
    elementSetName = config_json["coupling_params"]["element_set_to_couple"];
    std::pair<int, vector<double>> vertexInfo = getRelevantMaterialPoints(fem, elementSetName);
    vector<double> vertexPositions = vertexInfo.second;
    numberOfVerticies = vertexInfo.first;
    std::stringstream infoOutput;
    infoOutput << "Number of vertices variables: " << numberOfVerticies << "\n";
    infoOutput << "Dimensions: " << dimensions << "\n";
    infoOutput << "Coupling ElementSet with Name: " << elementSetName << "\n";
    double xMin, yMin, zMin = std::numeric_limits<double>::max();
    double xMax, yMax, zMax = std::numeric_limits<double>::min();
    for (int i = 0; i < vertexPositions.size(); i = i + 3) {
        double x = vertexPositions.at(i);
        double y = vertexPositions.at(i + 1);
        double z = vertexPositions.at(i + 2);
        if (x < xMin) {
            xMin = x;
        } else if (x > xMax) {
            xMax = x;
        }
        if (y < yMin) {
            yMin = y;
        } else if (y > yMax) {
            yMax = y;
        }
        if (z < zMin) {
            zMin = z;
        } else if (z > zMax) {
            zMax = z;
        }
    }
    infoOutput << "Mesh dimensions: [" << xMin << "," << yMin << ',' << zMin << "]x[" << xMax << "," << yMax << ","
               << zMax << "]\n";

    // precice write Setup
    std::string writeMeshName = config_json["coupling_params"]["write_mesh_name"];
    int writeMeshID = precice->getMeshID(writeMeshName);
    auto write_data_names = config_json["coupling_params"]["write_data_name"];
    infoOutput << "Number of write variables: " << write_data_names.size() << "\n";
    for (int i = 0; i < write_data_names.size(); i++) {
        std::string type = write_data_names[i]["type"];
        std::string mappingName = write_data_names[i]["mapping_name"];
        std::string febioClassName = write_data_names[i]["febio_class_name"];
        std::string febioVariableName = write_data_names[i]["name"];
        std::string preciceType = write_data_names[i]["precice_type"];
        std::string febioType = write_data_names[i]["febio_type"];
        infoOutput << mappingName.c_str() << "\n";
        int dataID = precice->getDataID(mappingName, writeMeshID);
        int vectorSize = 0;
        if (preciceType == "vector") {
            vectorSize = numberOfVerticies;
        } else {
            vectorSize = numberOfVerticies * dimensions;
        }

        FebioConfigEntry configEntry = {type, febioClassName, febioVariableName, preciceType, febioType, dataID,
                                        vector<double>(vectorSize, 0)};
        preciceWriteData.insert(std::pair<std::string, FebioConfigEntry>(mappingName, configEntry));
    }

    // preCICE read Setup
    std::string readMeshName = config_json["coupling_params"]["read_mesh_name"];
    int readMeshID = precice->getMeshID(readMeshName);
    auto read_data_names = config_json["coupling_params"]["read_data_name"];
    infoOutput << "Number of read variables: " << read_data_names.size() << "\n";
    for (int i = 0; i < read_data_names.size(); i++) {
        std::string type = read_data_names[i]["type"];
        std::string mappingName = read_data_names[i]["mapping_name"];
        std::string febioClassName = read_data_names[i]["febio_class_name"];
        std::string febioVariableName = read_data_names[i]["name"];
        std::string preciceType = read_data_names[i]["precice_type"];
        std::string febioType = read_data_names[i]["febio_type"];
        infoOutput << mappingName.c_str() << "\n";
        int dataID = precice->getDataID(mappingName, readMeshID);
        int vectorSize = 0;
        if (preciceType == "scalar") {
            vectorSize = numberOfVerticies;
        } else {
            vectorSize = numberOfVerticies * dimensions;
        }
        FebioConfigEntry configEntry = {type, febioClassName, febioVariableName, preciceType, febioType, dataID,
                                        vector<double>(vectorSize, 0)};
        preciceReadData.insert(std::pair<std::string, FebioConfigEntry>(mappingName, configEntry));
    }
    feLogInfo(infoOutput.str().c_str());

    // Setup vertices
    vertexIDs.resize(numberOfVerticies);
    precice->setMeshVertices(writeMeshID, numberOfVerticies, vertexPositions.data(), vertexIDs.data());


    precice_dt = precice->initialize();
    precice->initializeData();
    feLogInfo("Finished precice init");
}


bool PreciceCallback::Execute(FEModel &fem, int nreason) {
    if (nreason == CB_INIT) {
        Init(&fem);
    } else if (nreason == CB_UPDATE_TIME) {
        if (precice->isActionRequired(cowic)) {
            feLogInfo("CB_UPDATE_TIME - Saving Checkpoint\n");
            // Save
            // this uses dmp.open(true,true) which leads to the time controller not beeing serialized
            // Also setting dmp.open(true,false) leads to segfault dont know why yet
            // Switch time controller
            if (checkpointTimeStepController) { // Free earlier checkpointTimeStepContoller instance
                delete checkpointTimeStepController;
            }
            checkpointTimeStepController = fem.GetCurrentStep()->m_timeController;
            FETimeStepController *newTimeController = new FETimeStepController(&fem);
            newTimeController->SetAnalysis(fem.GetCurrentStep());
            newTimeController->CopyFrom(checkpointTimeStepController);
            fem.GetCurrentStep()->m_timeController = newTimeController;
            checkpoint_time = fem.GetTime().currentTime;
            dmp.clear();
            fem.Serialize(dmp);
            precice->markActionFulfilled(cowic);
        }
        dt = min(precice_dt, fem.GetCurrentStep()->m_dt);
        feLogInfo("Current Simulation Time %f\n", fem.GetTime().currentTime);
        feLogInfo("Timestep %f\n", dt);
        fem.GetCurrentStep()->m_dt = dt;
    } else if (nreason == CB_MAJOR_ITERS) {
        if (!precice->isCouplingOngoing()) {
            return false;
        } else {
            ReadData(&fem);
            WriteData(&fem);
            precice_dt = precice->advance(dt);
            if (precice->isActionRequired(coric)) {
                feLogInfo("CB_MAJOR_ITERS - Restoring Checkpoint\n");
                // Restore
                // taken from FEAnalysis.cpp Line 475 ff
                // restore the previous state
                dmp.Open(false, true); // This does not restore the time controller only if bshallow is false
                fem.Serialize(dmp);
                FETimeStepController *newTimeController = new FETimeStepController(&fem);
                newTimeController->SetAnalysis(fem.GetCurrentStep());
                newTimeController->CopyFrom(checkpointTimeStepController);
                fem.GetCurrentStep()->m_timeController = newTimeController;
                fem.GetTime().currentTime = checkpoint_time;
                fem.GetCurrentStep()->m_ntimesteps--; // Decrease number of steps because it gets increased right after this
                precice->markActionFulfilled(coric);
            }
        }
    } else if (nreason == CB_SOLVED) {
        precice->finalize();
        delete precice;
    }
    return true;
}

void PreciceCallback::ReadData(FEModel *fem) {
    if (precice->isReadDataAvailable()) {
        // Reading data from precice
        for (auto &[dataName, idVectorPair]: preciceReadData) {
            if (idVectorPair.preciceType == "scalar") {
                precice->readBlockScalarData(idVectorPair.dataID, numberOfVerticies, vertexIDs.data(),
                                             idVectorPair.data.data());
            } else if (idVectorPair.preciceType == "vector") {
                precice->readBlockVectorData(idVectorPair.dataID, numberOfVerticies, vertexIDs.data(),
                                             idVectorPair.data.data());
            } else {
                feLogError("Unknown preciceType");
                throw FEException("Unknown preciceType");
            }
        }

        int counter = 0;

        FEElementSet *elementSet = fem->GetMesh().FindElementSet(elementSetName);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);

                for (auto &[key, entry]: preciceReadData) {
                    std::variant<double, vector<double>> data;
                    if (entry.preciceType == "scalar") {
                        data = preciceReadData.at(key).data.at(counter);
                    } else if (entry.preciceType == "vector") {
                        vector<double> tmpData;
                        for (int k = 0; k <
                                        dimensions; k++) { // This should generate dimension matching vector with data from pecice
                            int index = counter * dimensions + k;
                            tmpData.push_back(preciceReadData.at(key).data.at(index));
                        }
                        data = tmpData;
                    }
                    insertData(materialPoint, entry.febioClassName, entry.febioVariableName,
                               entry.febioType, entry.preciceType, entry.type,
                               data);
                }
                counter++;
            }
        }
    }
}


void PreciceCallback::WriteData(FEModel *fem) {
    // Write data to precice
    if (precice->isWriteDataRequired(dt)) {
        FEMesh &femMesh = fem->GetMesh();
        // Updating Data that is send to the micromanager
        int counter = 0;
        FEElementSet *elementSet = fem->GetMesh().FindElementSet(elementSetName);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
                if (!materialPoint) {
                    continue;
                }

                for (auto &[key, entry]: preciceWriteData) {
                    auto value = extractData(materialPoint, entry.febioClassName, entry.febioVariableName,
                                             entry.febioType, entry.preciceType, entry.type);
                    if (holds_alternative<double>(value)) {
                        entry.data[counter] = std::get<double>(value);
                    } else if (holds_alternative<vector<double>>(value)) {
                        using std::begin, std::end;
                        vector<double> returnValue = std::get<vector<double>>(value);
                        entry.data.insert(end(entry.data), begin(returnValue), end(returnValue));
                    }
                }
                counter++;
            }
        }
        // Writing data to precice
        for (auto &[dataName, idVectorPair]: preciceWriteData) {
            if (idVectorPair.preciceType == "scalar") {
                precice->writeBlockScalarData(idVectorPair.dataID, numberOfVerticies, vertexIDs.data(),
                                              idVectorPair.data.data());
            } else if (idVectorPair.preciceType == "vector") {
                precice->writeBlockVectorData(idVectorPair.dataID, numberOfVerticies, vertexIDs.data(),
                                              idVectorPair.data.data());

            } else {
                feLogError("Unknown preciceType");
                throw FEException("Unknown preciceType");
            }
        }
    }
}

struct InsertDataFrompreCICEToFEBioWrapper {
public:
    FEMaterialPoint *materialPoint;
    const std::string &variableName;
    const std::string &preciceVariableType;
    const std::string &febioVariableType;
    const std::string &type;
    const std::string &className;
    const std::variant<double, vector<double>> &value;

    template<typename T>
    void insert(FEModel *model) {
        if (!materialPoint) {
            throw FEException("Material point is nullptr");
        }
        std::variant<vec3d, vec2d, double> valueExtracted;
        if (std::holds_alternative<double>(value)) {
            if (febioVariableType == "double") {
                valueExtracted = std::get<double>(value);
            } else if (febioVariableType == "vec2d") {
                vec2d valueExtractedTmp = {std::get<double>(value), 0};
                valueExtracted = valueExtractedTmp;
            } else if (febioVariableType == "vec3d") {
                vec3d valueExtractedTmp = {std::get<double>(value), 0, 0};
                valueExtracted = valueExtractedTmp;
            } else {
                model->Logf(2, "Variable Conversion is not implemented");
                throw FEException("Variable Conversion is not implemented");
            }
        } else if (std::holds_alternative<vector<double>>(value)) {
            if (febioVariableType == "double") {
                valueExtracted = std::get<vector<double>>(value)[0];
            } else if (febioVariableType == "vec2d") {
                vector<double> vectorExtracted = std::get<vector<double>>(value);
                vec2d valueExtractedTmp = {vectorExtracted[0], vectorExtracted[1]};
                valueExtracted = valueExtractedTmp;
            } else if (febioVariableType == "vec3d") {
                vector<double> vectorExtracted = std::get<vector<double>>(value);
                vec3d valueExtractedTmp = {vectorExtracted[0], vectorExtracted[1], vectorExtracted[2]};
                valueExtracted = valueExtractedTmp;
            } else {
                model->Logf(2, "Variable Conversion is not implemented");
                throw FEException("Variable Conversion is not implemented");
            }

        }

        T *castPoint = dynamic_cast<T *>(materialPoint);
        if (type == "variable") {
            property prop = type::get(castPoint).get_property(variableName);
            if (prop.is_valid()) {
                prop.set_value(castPoint, valueExtracted);
            } else {
                model->Logf(2, "Could not find class Property, is it register with rttr?");
            }
        } else if (type == "function") {
            method meth = type::get(castPoint).get_method(variableName);
            if (meth.is_valid()) {
                if (holds_alternative<double>(valueExtracted)) {
                    meth.invoke(castPoint, get<double>(valueExtracted));
                } else {
                    model->Logf(2, "Not implemented fix this");
                }
            } else {
                model->Logf(2, (std::string("FEBio Method not found: ") + variableName +
                                std::string(" on variable type: ") + febioVariableType + std::string(" type: ") +
                                type + std::string(" className: ") + className).c_str());
                throw FEException("FEBio Method not found");
            }
        }
    }
};

void PreciceCallback::insertData(FEMaterialPoint *materialPoint, const std::string &className,
                                 const std::string &variableName,
                                 const std::string &febioVariableType, const std::string &preciceVariableType,
                                 const std::string &type,
                                 const std::variant<double, vector<double>> &value) {

    InsertDataFrompreCICEToFEBioWrapper wrapper = {materialPoint, variableName, preciceVariableType, febioVariableType,
                                                   type,
                                                   className, value};
    if (className == "FEBiphasicMaterialPointTPM") {
        wrapper.insert<FEBiphasicMaterialPointTPM>(GetFEModel());
    } else if (className == "FESolutesMaterialPointTPM") {
        wrapper.insert<FESolutesMaterialPointTPM>(GetFEModel());
    } else if (className == "FEFluidMaterialPoint") {
        wrapper.insert<FEFluidMaterialPoint>(GetFEModel());
    } else {
        feLogError((std::string("Unknown Class Name ") + className).c_str());
    }


}

struct ExtractDataFromFEBioWrapper {
public:
    FEMaterialPoint *materialPoint;
    const std::string &variableName;
    const std::string &preciceVariableType;
    const std::string &febioVariableType;
    const std::string &type;
    const std::string &className;

    template<typename T>
    std::variant<double, vector<double>> extract(FEModel *model) {
        rttr::variant var_prop;
        if (!materialPoint) {
            throw FEException("Material point is nullptr");
        }
        const T &pt = *materialPoint->ExtractData<T>();
        if (type == "variable") {
            property prop = type::get(pt).get_property(variableName);
            if (prop.is_valid()) {
                var_prop = prop.get_value(pt);
            } else {
                model->Logf(2, "Could not find class Property, is it register with rttr?");
            }
        } else if (type == "function") {
            method meth = type::get(pt).get_method(variableName);
            if (meth) {
                var_prop = meth.invoke(pt);
            } else {
                model->Logf(2, (std::string("FEBio Method not found: ") + variableName +
                                std::string(" on variable type: ") + febioVariableType + std::string(" type: ") +
                                type + std::string(" className: ") + className).c_str());
                throw FEException("FEBio Method not found");
            }
        } else {
            model->Logf(2, (std::string("Unknown type ") + type).c_str());
            throw FEException("Unknown type either function or variable");
        }
        if (preciceVariableType == "scalar") {
            if (febioVariableType == "double") {
                double value = var_prop.get_value<double>();
                return value;
            } else if (febioVariableType == "vector<double>") {
                vector<double> value = var_prop.get_value<vector<double>>();
                return value.at(0); // TODO this has to be documented !!!!
            } else {
                model->Logf(2, (std::string("FEBio Datatype not implemnted/unknown: ") + febioVariableType).c_str());
                throw FEException("FEBio Datatype not implemented/unknown");
            }

        } else if (preciceVariableType == "vector") {
            if (febioVariableType == "double") {
                vector<double> value = {var_prop.get_value<double>()};
                return value;
            } else if (febioVariableType == "vector<double>") {
                vector<double> value = var_prop.get_value<vector<double>>();
                return value;
            } else if (febioVariableType == "vec2d") {
                vec2d tempValue = var_prop.get_value<vec2d>();
                vector<double> value = {tempValue.x(), tempValue.y()};
                return value;
            } else if (febioVariableType == "vec3d") {
                vec3d tempValue = var_prop.get_value<vec3d>();
                vector<double> value = {tempValue.x, tempValue.y, tempValue.z};
                return value;
            } else {
                model->Logf(2, (std::string("FEBio Datatype not implemnted/unknown: ") + febioVariableType).c_str());
                throw FEException("FEBio Datatype not implemented/unknown");
            }
        } else {
            model->Logf(2, (std::string("preCICEVariableType unknown: ") + preciceVariableType).c_str());
            throw FEException("preciceVariableType Unknown");
        }
    }
};

std::variant<double, vector<double>>
PreciceCallback::extractData(FEMaterialPoint *materialPoint, const std::string &className,
                             const std::string &variableName,
                             const std::string &febioVariableType, const std::string &preciceVariableType,
                             const std::string &type) {
    ExtractDataFromFEBioWrapper wrapper = {materialPoint, variableName, preciceVariableType, febioVariableType, type,
                                           className};
    if (className == "FEBiphasicMaterialPointTPM") {
        return wrapper.extract<FEBiphasicMaterialPointTPM>(GetFEModel());
    } else if (className == "FESolutesMaterialPointTPM") {
        return wrapper.extract<FESolutesMaterialPointTPM>(GetFEModel());
    } else if (className == "FEFluidMaterialPoint") {
        return wrapper.extract<FEFluidMaterialPoint>(GetFEModel());
    } else {
        feLogError((std::string("Unknown Class Name ") + className).c_str());
    }

}
