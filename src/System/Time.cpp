#include "Time.h"

void Time::start()
{
    start_ = std::chrono::steady_clock::now();
}

void Time::stop()
{
    end_ = std::chrono::steady_clock::now();
}

double Time::elapsedSeconds() const
{
    return (end_ - start_).count()/1e9;
}
