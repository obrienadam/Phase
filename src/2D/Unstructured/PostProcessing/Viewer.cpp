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

    split(integerFields_, integerFields, is_any_of(", "), token_compress_on);
    split(scalarFields_, scalarFields, is_any_of(", "), token_compress_on);
    split(vectorFields_, vectorFields, is_any_of(", "), token_compress_on);
}
