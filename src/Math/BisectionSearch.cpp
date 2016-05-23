#include <math.h>

#include "BisectionSearch.h"

std::pair<Scalar, bool> bisectionSearch(const std::function<Scalar(Scalar)>& f, Scalar a, Scalar b, Scalar toler, int maxIters)
{
    Scalar x;

    for(int i = 0; i < maxIters; ++i)
    {
        x = (a + b)/2.;
        Scalar result = f(x);

        if(fabs(result) <= toler)
            return std::make_pair(x, true);
        else if(isnan(result))
            break;
        else if (result > 0.)
            b = x;
        else
            a = x;
    }

    return std::make_pair(x, false);
}
