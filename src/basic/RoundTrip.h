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
struct round_trip
{
    explicit round_trip(double v) : value(v) {}
    double value;
};

inline std::ostream & operator<<(std::ostream & out, const round_trip & r)
{
    char buf[32];
    for (int digits = 15; digits <= 17; digits++)
    {
        snprintf(buf, sizeof(buf), "%.*g", digits, r.value);
        if (strtod(buf, NULL) == r.value) break;
    }
    return out << buf;
}

/// real numbers are written to round trip; other values as they are
inline round_trip round_trip_value(double v)
{
    return round_trip(v);
}
inline round_trip round_trip_value(float v)
{
    return round_trip(v);
}
template <typename T> inline const T & round_trip_value(const T & v)
{
    return v;
}

} // namespace cusg
