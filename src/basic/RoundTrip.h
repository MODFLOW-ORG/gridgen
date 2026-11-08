/*! \file RoundTrip.h
\brief write a real number with the digits needed to read it back exactly
*/
#pragma once

#include <cstdio>
#include <cstdlib>
#include <ostream>

namespace cusg
{

/// a real number that is written with the fewest significant digits, up to
/// 17, that read back to the same value
struct RoundTrip
{
    explicit RoundTrip(double v) : value(v) {}
    double value;
};

inline std::ostream & operator<<(std::ostream & out, const RoundTrip & r)
{
    char buf[32];
    for (int digits = 15; digits <= 17; digits++)
    {
        snprintf(buf, sizeof(buf), "%.*g", digits, r.value);
        if (strtod(buf, nullptr) == r.value) break;
    }
    return out << buf;
}

/// real numbers are written to round trip; other values as they are
inline RoundTrip round_trip_value(double v)
{
    return RoundTrip(v);
}
inline RoundTrip round_trip_value(float v)
{
    return RoundTrip(v);
}
template <typename T> inline const T & round_trip_value(const T & v)
{
    return v;
}

} // namespace cusg
