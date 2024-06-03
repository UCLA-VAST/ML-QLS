#ifndef BITWUZLATERMINATOR_HPP
#define BITWUZLATERMINATOR_HPP

#include <bitwuzla/cpp/bitwuzla.h>

#include <chrono>

using namespace bitwuzla;
using namespace std::chrono;


class RuntimeTerminator : public Terminator
{
 public:
  RuntimeTerminator(uint32_t time_limit_ms)
      : Terminator(),
        time_limit_ms(time_limit_ms),
        start(high_resolution_clock::now())
  {
  }
  bool terminate() override
  {
    if (duration_cast<milliseconds>(high_resolution_clock::now() - start)
            .count()
        >= time_limit_ms)
    {
      return true;
    }
    return false;
  }
  uint32_t time_limit_ms = 0;
  high_resolution_clock::time_point start;
};

#endif