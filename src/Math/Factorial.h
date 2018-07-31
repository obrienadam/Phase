#ifndef PHASE_FACTORIAL_H
#define PHASE_FACTORIAL_H

template<typename T>
T factorial(T n)
{
    T r = 1;

    for(T i = 2; i <= n; ++i)
        r *= i;

    return r;
}

#endif
