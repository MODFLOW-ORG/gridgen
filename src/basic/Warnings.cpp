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

void warn_features_outside(int noutside, int nfeatures, const std::string & shapefile,
                           const std::string & block, int layer, const std::string & effect)
{
    if (noutside == 0) return;
    std::string source = " in shapefile " + shapefile + " (" + block + ") ";
    std::string subject;
    if (nfeatures == 1)
        subject = "The feature" + source + "does not intersect";
    else if (noutside == nfeatures)
        subject = "None of the " + std::to_string(nfeatures) + " features" + source + "intersect";
    else
        subject = std::to_string(noutside) + " of the " + std::to_string(nfeatures) + " features" +
                  source + "do not intersect";
    //the effect refers to the features as "them"
    std::string consequence = effect;
    if (noutside == 1)
    {
        size_t pos = consequence.find("them");
        if (pos != std::string::npos) consequence.replace(pos, 4, "it");
    }
    warn(subject + " layer " + std::to_string(layer + 1) + " of the grid; " + consequence + ".");
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
