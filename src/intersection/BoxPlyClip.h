#pragma once
#include "intersection.h"
#include "polygon.h"
#include <list>
#include <vector>

using namespace std;

#define EPSILON 0.0001
#define SMALLCOORDINATE 1e-8
#define SMALLAREA 1e-7

namespace cusg
{
class BoxPlyClip;

struct VNode
{
    double x, y;
    struct VNode * next;
    struct VNode * prev;
    bool isInDirect;
    bool isOnline;
    bool isVisited;
    bool isScanPoint;
    VNode() { init(); }
    VNode(double x_coordinate, double y_coordinate)
    {
        x = x_coordinate;
        y = y_coordinate;
        init();
    }
    ~VNode()
    {
        next = nullptr;
        prev = nullptr;
    }
    void init()
    {
        next = nullptr;
        prev = nullptr;
        isInDirect = false;
        isOnline = false;
        isVisited = false;
        isScanPoint = false;
    }
    void resetStatus()
    {
        isInDirect = false;
        isOnline = false;
        isVisited = false;
        isScanPoint = false;
    }
    void setDirection(double line_coordinate, bool isX, bool isGreater)
    {
        if (isX)
        {
            if (fabs(x - line_coordinate) < SMALLNUMBER)
            {
                isInDirect = true;
                return;
            }
            bool big = (x >= line_coordinate);
            isInDirect = (big == isGreater);
        }
        else
        {
            if (fabs(y - line_coordinate) < SMALLNUMBER)
            {
                isInDirect = true;
                return;
            }
            bool big = (y >= line_coordinate);
            isInDirect = (big == isGreater);
        }
    }
    bool Online(double line_coordinate, bool isX)
    {
        if (isX && (fabs(x - line_coordinate) < SMALLNUMBER))
            return true;
        else if ((!isX) && (fabs(y - line_coordinate) < SMALLNUMBER))
            return true;
        return false;
    }
};

class BoxPlyClip
{
public:
    BoxPlyClip(void);
    ~BoxPlyClip(void);

    void destroy();
    void destroy(VNode * head);
    void resetStatus();

    //add a vertex, which this clip then owns
    void add(VNode * v);
    //insert a vertex
    void insert(VNode * before, VNode * after, VNode * v);
    //next
    VNode * next(VNode * v);

    //clip using scanning line
    void clip(double line_coordinate, bool isX, bool isBigger, list<BoxPlyClip *> & clips);
    //clip using a box
    void clip(double xmin, double xmax, double ymin, double ymax, list<BoxPlyClip *> & clips);

    //set intersection based on the scanning line
    int setIntersection(VNode * start_vertex, double line_coordinate, bool isX, bool isBigger,
                        vector<VNode *> & line_vertices);

    //collect
    void collectPolygonsInInterval(bool isX, bool isBigger, vector<VNode *> & line_vertices,
                                   list<BoxPlyClip *> & clips, int from, int to, bool forward);

    //collect polygon results
    void collectPolygons(double line_coordinate, bool isX, bool isBigger,
                         vector<int> & line_vertex_counts, vector<VNode *> & line_vertices,
                         list<BoxPlyClip *> & clips);
    //collect one ply
    //BoxPlyClip* collectOnePly(vector<VNode*>& line_vertices, int start_index, int end_index, bool isNextOrder);
    void collectOnePly(vector<VNode *> & line_vertices, int start_index, int end_index,
                       bool isNextOrder, list<BoxPlyClip *> & clip_results, bool isX,
                       bool isBigger);

    void collectSecondPart(vector<VNode *> & line_vertices, int start_index, int end_index,
                           VNode * start_node, VNode * end_node, bool forward, int vertex_count,
                           list<BoxPlyClip *> & clip_results, bool isX, bool isBigger,
                           vector<VNode *> & secondHalf);

    void constructClipPly(vector<VNode *> & firstHalf, vector<VNode *> & secondHalf,
                          list<BoxPlyClip *> & clip_results);

    VNode * collectAlongWay(vector<VNode *> & vertices, VNode * start, bool forward);
    //unprocessed
    //bool unprocessed();
    VNode * copy(const VNode * v);

    //this vector stores the array of holes
    VNode * first;      //mark the outboundary
    list<VNode *> plys; //make sure the order, out boundary is the first one

    //new added vertices
    list<VNode *> added_vertices;
    //influenced vertices
    list<VNode *> changed_vertices;
};
//convert the ply to BoxPlyClip
BoxPlyClip * convertPlyToClip(const c_ply & ply);

//convert BoxPlyClip to c_ply
//void convertClipToPly(BoxPlyClip* tcp, c_ply& ply);

void convertClipToPly(VNode * first_vertex, c_ply & ply);

BoxPlyClip * convertCPolygonToClip(const GIS_polygon & polygon);

BoxPlyClip * convertSimplePolygonToClip(const c_polygon & polygon);

void convertClipToPolygon(BoxPlyClip * ply_clip, GIS_polygon * polygon);

VNode * copyOnePly(VNode * head, bool bSetVisited = false);

double getRotAngle(Vector2d & v1, Vector2d & v2);

//make sure the min and max are in correct order
double intersectTwoBox(Point2d lower_left, Point2d lower_right, Point2d upper_right,
                       Point2d upper_left, double xmin2, double xmax2, double ymin2, double ymax2);

}