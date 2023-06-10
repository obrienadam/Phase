#include <sstream>

#include "Timer.h"

void Timer::start() { start_ = std::chrono::steady_clock::now(); }

void Timer::stop() { end_ = std::chrono::steady_clock::now(); }

void Timer::reset() { start_ = end_ = std::chrono::steady_clock::now(); }

double Timer::elapsedSeconds() const {
  return std::chrono::duration_cast<Seconds>(end_ - start_).count();
}

double Timer::elapsedSeconds(const Communicator &comm) const {
  return comm.broadcast(comm.mainProcNo(), elapsedSeconds());
}

std::string Timer::elapsedTime() const {
  auto dur(end_ - start_);
  std::ostringstream sout;

  Hours hours = std::chrono::duration_cast<Hours>(dur);
  dur -= hours;

  Minutes minutes = std::chrono::duration_cast<Minutes>(dur);
  dur -= minutes;

  Seconds seconds = std::chrono::duration_cast<Seconds>(dur);

  sout << hours.count() << ":" << minutes.count() << ":" << seconds.count();

  return sout.str();
}

std::string Timer::elapsedCpuTime(const Communicator &comm) const {
  std::chrono::duration<double, std::nano> dur(
      comm.sum((end_ - start_).count()));
  std::ostringstream sout;

  Hours hours = std::chrono::duration_cast<Hours>(dur);
  dur -= hours;

  Minutes minutes = std::chrono::duration_cast<Minutes>(dur);
  dur -= minutes;

  Seconds seconds = std::chrono::duration_cast<Seconds>(dur);

  sout << hours.count() << ":" << minutes.count() << ":" << seconds.count();

  return sout.str();
}
