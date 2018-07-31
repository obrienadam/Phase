#include "System/Exception.h"

#include "SparseMatrixSolver.h"

void SparseMatrixSolver::set(const std::vector<std::tuple<Index, Index, Scalar>> &entries)
{
    CoefficientList coeffs;
    coeffs.reserve(entries.size() / 5);

    for (const auto &entry: entries)
    {
        Index row = std::get<0>(entry);
        Index col = std::get<1>(entry);
        Scalar val = std::get<2>(entry);

        if (row >= coeffs.size())
            coeffs.resize(row + 1);

        bool newEntry = true;
        for(Entry& entry: coeffs[row])
            if(entry.first == col)
            {
                entry.second += val;
                newEntry = false;
                break;
            }

        if(newEntry)
            coeffs[row].push_back(Entry(col, val));
    }

    set(coeffs);
}

Scalar SparseMatrixSolver::solve(const Vector &x0)
{
    setGuess(x0);
    return solve();
}

Scalar SparseMatrixSolver::solveLeastSquares()
{
    throw Exception("SparseMatrixSolver", "solveLeastSquares", "least squares solver is not available for this sparse matrix solver type.");
}

void SparseMatrixSolver::printStatus(const std::string &msg) const
{
    printf("%s iterations = %d, error = %lf.\n", msg.c_str(), nIters(), error());
}
