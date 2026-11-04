#include "Warnings.h"

#include <iostream>

namespace cusg
{

namespace
{
std::vector<std::string> recorded;
}

void warn(const std::string & message)
{
    std::cerr << "! Warning: " << message << std::endl;
    recorded.push_back(message);
}

void warn_features_outside(int count, int total, const std::string & shapefile,
                           const std::string & block, int layer, const std::string & effect)
{
    if (count == 0) return;
    std::string where = " in shapefile " + shapefile + " (" + block + ") ";
    std::string what;
    if (total == 1)
        what = "The feature" + where + "does not intersect";
    else if (count == total)
        what = "None of the " + std::to_string(total) + " features" + where + "intersect";
    else
        what = std::to_string(count) + " of the " + std::to_string(total) + " features" + where +
               "do not intersect";
    //the effect refers to the features as "them"
    std::string said = effect;
    if (count == 1)
    {
        size_t pos = said.find("them");
        if (pos != std::string::npos) said.replace(pos, 4, "it");
    }
    warn(what + " layer " + std::to_string(layer + 1) + " of the grid; " + said + ".");
}

const std::vector<std::string> & get_warnings()
{
    return recorded;
}

void clear_warnings()
{
    recorded.clear();
}

} // namespace cusg
