#include <boost/filesystem.hpp>

#include "Viewer.h"

Viewer::Viewer(const Input &input, const Communicator &comm, const Solver &solver)
    :
      solver_(solver)
{
    using namespace std;
    using namespace boost;

    filename_ = input.caseInput().get<std::string>("CaseName");

    string vectorFields = input.caseInput().get<string>("Viewer.vectorFields");
    string scalarFields = input.caseInput().get<string>("Viewer.scalarFields");

    vector<string> vectorFieldNames, scalarFieldNames;
    split(vectorFieldNames, vectorFields, is_any_of(", "), token_compress_on);
    split(scalarFieldNames, scalarFields, is_any_of(", "), token_compress_on);

    for(const auto& field: solver_.scalarFields())
    {
        if(std::find(scalarFieldNames.begin(), scalarFieldNames.end(), field.first) != scalarFieldNames.end())
            scalarFields_.push_back(Ref<const ScalarFiniteVolumeField>(field.second));
    }

    for(const auto& field: solver_.vectorFields())
    {
        if(std::find(vectorFieldNames.begin(), vectorFieldNames.end(), field.first) != vectorFieldNames.end())
            vectorFields_.push_back(Ref<const VectorFiniteVolumeField>(field.second));
    }
}
