//
// Created by fritz on 2/25/22.
//

#include <map>
#include <FECore/sdk.h>
#include <FECore/FECallBack.h>
#include <FECore/Callback.h>
#include <FECore/FENode.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEAnalysis.h>
#include <precice/SolverInterface.hpp>
#include <FECore/DumpMemStream.h>
#include <rttr/variant.h>
#include <variant>

class PreciceCallback : public FECallBack {
public:

    PreciceCallback(FEModel *pfem) : FECallBack(pfem, CB_INIT | CB_UPDATE_TIME | CB_MAJOR_ITERS), dmp(*pfem) {};

    void ReadData(FEModel *fem);

    void WriteData(FEModel *fem);

    void Init(FEModel *fem);

    bool Execute(FEModel &fem, int nreason);

    std::variant<double, vector<double>>
    extractData(FEMaterialPoint *materialPoint, const std::string &className,
                const std::string &variableName, const std::string &febioVariableType,
                const std::string &preciceVariableType, const std::string &type);

    void insertData(FEMaterialPoint *materialPoint, const std::string &className,
                    const std::string &variableName, const std::string &febioVariableType,
                    const std::string &preciceVariableType, const std::string &type,
                    const std::variant<double, vector<double>> &value);

    std::pair<int, vector<double>> getRelevantMaterialPoints(FEModel* fem, const std::string &elementName);

protected:
    precice::SolverInterface *precice = NULL;
    const std::string &coric = precice::constants::actionReadIterationCheckpoint();
    const std::string &cowic = precice::constants::actionWriteIterationCheckpoint();
    FEAnalysis *checkPointStep;
    double dimensions; // preCICE dimensions
    double dt; // solver timestep size
    double precice_dt; // maximum precice timestep size
    int numberOfVerticies; // number of vertices at wet surface
    std::vector<int> vertexIDs;
    std::string elementSetName;

    struct FebioConfigEntry {
        std::string type;
        std::string febioClassName;
        std::string febioVariableName;
        std::string preciceType;
        std::string febioType;
        int dataID;
        vector<double> data;
    };
    // {precice_variable_name, {febioClassName, febioVariableName, preciceType, febioType, dataID, vectorIDs}}
    std::map<std::string, FebioConfigEntry> preciceReadData;
    std::map<std::string, FebioConfigEntry> preciceWriteData;

    DumpMemStream dmp;
    double checkpoint_time = 0;
    FETimeStepController *checkpointTimeStepController = nullptr;

};


