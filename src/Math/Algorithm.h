#ifndef PHASE_ALGORITHM_H
#define PHASE_ALGORITHM_H

#include <algorithm>

template<class T>
bool inRange(const T& val, const T& low, const T& high)
{
    return (val >= low) && (val <= high);
}

template<class T>
constexpr const T& clamp(const T& val, const T& low, const T& high)
{
    return std::max(std::min(val, high), low);
}

template<class T, class TFunc>
T bisectionSearch(std::pair<T, T> bounds, const TFunc& func, const T& tolerance, size_t maxIters)
{
    T x = (bounds.first + bounds.second) / 2.;
    T y = func(x);

    if(func(bounds.second) < func(bounds.first))
        std::swap(bounds.first, bounds.second);

    for(size_t iter = 0; iter < maxIters && std::abs(y) > tolerance; ++iter)
    {
        if(y > 0)
            bounds.second = x;
        else
            bounds.first = x;

        x = (bounds.first + bounds.second) / 2.;
        y = func(x);
    }

    return x;
}

template<class T, class TFunc>
T secantSearch(std::pair<T, T> x, const TFunc& func, const T& tolerance, size_t maxIters)
{
    std::pair<T, T> f = std::make_pair(
            func(x.first),
            func(x.second)
    );

    for(size_t iter = 0; iter < maxIters && std::abs(f.first) > tolerance; ++iter)
    {
        x = std::make_pair((x.second*f.first - x.first*f.second) / (f.first - f.second), x.first);
        f = std::make_pair(func(x.first), f.first);
    }

    return x.second;
}

#endif
