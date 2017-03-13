#ifndef BISECTION_SEARCH_H
#define BISECTION_SEARCH_H

#include <functional>

#include "Types.h"

std::pair<Scalar, bool>
bisectionSearch(const std::function<Scalar(Scalar)> &f, Scalar a, Scalar b, Scalar toler = 1e-12, int maxIters = 100);

#endif
