#ifndef TIME_H
#define TIME_H

#include <chrono>

#include "Communicator.h"

class Time
{
public:

    typedef std::chrono::hours Hours;
    typedef std::chrono::minutes Minutes;
    typedef std::chrono::duration<double, std::ratio<1>> Seconds;

    void start();
    void stop();
    void reset();

    double elapsedSeconds() const;
    double elapsedSeconds(const Communicator& comm) const;

    std::string elapsedTime() const;
    std::string elapsedCpuTime(const Communicator& comm) const;

private:

    std::chrono::time_point<std::chrono::steady_clock> start_, end_;

};

#endif
