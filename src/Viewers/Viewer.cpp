#include <boost/filesystem.hpp>

#include "Viewer.h"

Viewer::Viewer(const Input &input, const Solver &solver)
    :
      solver_(solver)
{
    using namespace std;
    using namespace boost;

    filename_ = input.caseInput().get<std::string>("CaseName");

    string integerFields = input.caseInput().get<string>("Viewer.integerFields", "");
    string scalarFields = input.caseInput().get<string>("Viewer.scalarFields", "");
    string vectorFields = input.caseInput().get<string>("Viewer.vectorFields", "");

    vector<string> integerFieldNames, scalarFieldNames, vectorFieldNames;
    split(integerFieldNames, integerFields, is_any_of(", "), token_compress_on);
    split(scalarFieldNames, scalarFields, is_any_of(", "), token_compress_on);
    split(vectorFieldNames, vectorFields, is_any_of(", "), token_compress_on);

    for(const auto& field: solver_.integerFields())
    {
        if(std::find(integerFieldNames.begin(), integerFieldNames.end(), field.first) != integerFieldNames.end())
            integerFields_.push_back(Ref<const FiniteVolumeField<int> >(*field.second));
    }

    for(const auto& field: solver_.scalarFields())
    {
        if(std::find(scalarFieldNames.begin(), scalarFieldNames.end(), field.first) != scalarFieldNames.end())
            scalarFields_.push_back(Ref<const ScalarFiniteVolumeField>(*field.second));
    }

    for(const auto& field: solver_.vectorFields())
    {
        if(std::find(vectorFieldNames.begin(), vectorFieldNames.end(), field.first) != vectorFieldNames.end())
            vectorFields_.push_back(Ref<const VectorFiniteVolumeField>(*field.second));
    }
}
