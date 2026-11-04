/*! \file Warnings.h
\brief warnings about input that is ignored or has no effect
*/
#pragma once

#include <string>
#include <vector>

namespace cusg
{

/// print a warning to stderr and record it
void warn(const std::string & message);

/// warn that some features of a shapefile intersect no cell, and say what
/// that means for the grid
void warn_features_outside(int count, int total, const std::string & shapefile,
                           const std::string & block, int layer, const std::string & effect);

/// the warnings recorded since the last clear_warnings()
const std::vector<std::string> & get_warnings();

/// forget the recorded warnings
void clear_warnings();

} // namespace cusg
