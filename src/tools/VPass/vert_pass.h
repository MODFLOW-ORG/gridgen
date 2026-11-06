
#pragma once

#include "Grid.h"
#include "ModflowGrid2D.h"
#include "ModflowGrid.h"
#include "QuadTree3D.h"
#include <list>

using namespace std;
using namespace cusg;

struct VPTPair
{
    VPTPair()
    {
        source = nullptr;
        destination = nullptr;
        area = 0.0;
    }
    VPTPair(Box * source_box, Box * destination_box, double pass_area)
    {
        source = source_box;
        destination = destination_box;
        area = pass_area;
    }
    VPTPair(const VPTPair & other)
    {
        source = other.source;
        destination = other.destination;
        area = other.area;
    }
    Box * source;
    Box * destination;
    double area;
};

struct VPassBox
{
    VPassBox()
    {
        source = nullptr;
        destination = nullptr;
        area = 0.0;
    }
    VPassBox(const VPassBox & other)
    {
        source = other.source;
        destination = other.destination;
        pass_through.clear();
        pass_through.insert(pass_through.end(), other.pass_through.begin(),
                            other.pass_through.end());
        area = other.area;
    }

    Box * source;
    Box * destination;
    list<Box *> pass_through;
    double area;
};

void getVPT(Box * b, vector<VPTPair> & pairs, bool topToDown);
