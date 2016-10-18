#ifndef TIME_H
#define TIME_H

#include<chrono>

class Time
{
public:

    void start();
    void stop();

    double elapsedSeconds() const;

private:

    std::chrono::time_point<std::chrono::steady_clock> start_, end_;

};

#endif
