#ifndef FLUIDSIM_TIMESTYPES_H
#define FLUIDSIM_TIMESTYPES_H

#include <chrono>

namespace timetypes{

    using ns = std::chrono::nanoseconds;
    using sclock = std::chrono::steady_clock;
    using tp = sclock::time_point;

}




#endif //FLUIDSIM_TIMESTYPES_H
