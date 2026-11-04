/*! \file GridgenError.h
\brief error that stops a gridgen action
*/
#pragma once

#include <stdexcept>
#include <string>

namespace cusg
{

/// an error in the input or the grid that stops the action; the caller
/// reports the message
class GridgenError : public std::runtime_error
{
public:
    explicit GridgenError(const std::string & message) : std::runtime_error(message) {}
};

} // namespace cusg
