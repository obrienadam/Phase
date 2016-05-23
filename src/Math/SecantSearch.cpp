#include <math.h>

#include "SecantSearch.h"
#include "BisectionSearch.h"

std::pair<Scalar, bool> secantSearch(const std::function<Scalar(Scalar)>& f, Scalar a, Scalar b, Scalar toler, int maxIters)
{
    Scalar x = (a + b)/2., x0 = bisectionSearch(f, a, b, 0., 1).first, x1;
    Scalar result = f(x), result0 = f(x0);

    for(int i = 0; i < maxIters; ++i)
    {
        x1 = x - result/((result - result0)/(x - x0));

        x0 = x;
        x = x1;
        result0 = result;
        result = f(x);

        if(fabs(result) < toler)
            return std::make_pair(x, true);
        else if(isnan(x))
            return std::make_pair(x0, false);
        else if(isnan(result))
            break;
    }

    return std::make_pair(x, false);
}
