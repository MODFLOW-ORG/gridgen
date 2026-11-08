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
class RingClipper;

struct VNode
{
    double x, y;
    struct VNode * next;
    struct VNode * prev;
    bool is_in_direct;
    bool is_online;
    bool is_visited;
    bool is_scan_point;
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
        is_in_direct = false;
        is_online = false;
        is_visited = false;
        is_scan_point = false;
    }
    void reset_status()
    {
        is_in_direct = false;
        is_online = false;
        is_visited = false;
        is_scan_point = false;
    }
    void set_direction(double line_coordinate, bool is_x, bool is_greater)
    {
        if (is_x)
        {
            if (fabs(x - line_coordinate) < SMALLNUMBER)
            {
                is_in_direct = true;
                return;
            }
            bool big = (x >= line_coordinate);
            is_in_direct = (big == is_greater);
        }
        else
        {
            if (fabs(y - line_coordinate) < SMALLNUMBER)
            {
                is_in_direct = true;
                return;
            }
            bool big = (y >= line_coordinate);
            is_in_direct = (big == is_greater);
        }
    }
    bool on_line(double line_coordinate, bool is_x)
    {
        if (is_x && (fabs(x - line_coordinate) < SMALLNUMBER))
            return true;
        else if ((!is_x) && (fabs(y - line_coordinate) < SMALLNUMBER))
            return true;
        return false;
    }
};

class RingClipper
{
public:
    RingClipper(void);
    ~RingClipper(void);

    void destroy();
    void destroy(VNode * head);
    void reset_status();

    //add a vertex, which this clip then owns
    void add(VNode * v);
    //insert a vertex
    void insert(VNode * before, VNode * after, VNode * v);
    //next
    VNode * next(VNode * v);

    //clip using scanning line
    void clip(double line_coordinate, bool is_x, bool is_bigger, list<RingClipper *> & clips);
    //clip using a box
    void clip(double xmin, double xmax, double ymin, double ymax, list<RingClipper *> & clips);

    //set intersection based on the scanning line
    int set_intersection(VNode * start_vertex, double line_coordinate, bool is_x, bool is_bigger,
                         vector<VNode *> & line_vertices);

    //collect
    void collect_polygons_in_interval(bool is_x, bool is_bigger, vector<VNode *> & line_vertices,
                                      list<RingClipper *> & clips, int from, int to, bool forward);

    //collect polygon results
    void collect_polygons(double line_coordinate, bool is_x, bool is_bigger,
                          vector<int> & line_vertex_counts, vector<VNode *> & line_vertices,
                          list<RingClipper *> & clips);
    //collect one ring
    //RingClipper* collect_ring(vector<VNode*>& line_vertices, int start_index, int end_index, bool is_next_order);
    void collect_ring(vector<VNode *> & line_vertices, int start_index, int end_index,
                      bool is_next_order, list<RingClipper *> & clip_results, bool is_x,
                      bool is_bigger);

    void collect_second_part(vector<VNode *> & line_vertices, int start_index, int end_index,
                             VNode * start_node, VNode * end_node, bool forward, int vertex_count,
                             list<RingClipper *> & clip_results, bool is_x, bool is_bigger,
                             vector<VNode *> & second_half);

    void construct_clip_rings(vector<VNode *> & first_half, vector<VNode *> & second_half,
                              list<RingClipper *> & clip_results);

    VNode * collect_along_way(vector<VNode *> & vertices, VNode * start, bool forward);
    //unprocessed
    //bool unprocessed();
    VNode * copy(const VNode * v);

    //this vector stores the array of holes
    VNode * first;       //mark the outboundary
    list<VNode *> rings; //make sure the order, out boundary is the first one

    //new added vertices
    list<VNode *> added_vertices;
    //influenced vertices
    list<VNode *> changed_vertices;
};
//convert the ring to RingClipper
RingClipper * convert_ring_to_clip(const Ring & ring);

//convert RingClipper to Ring
//void convert_clip_to_ring(RingClipper* tcp, Ring& ring);

void convert_clip_to_ring(VNode * first_vertex, Ring & ring);

RingClipper * convert_polygon_to_clip(const GisPolygon & polygon);

RingClipper * convert_simple_polygon_to_clip(const Polygon & polygon);

void convert_clip_to_polygon(RingClipper * ring_clip, GisPolygon * polygon);

VNode * copy_ring(VNode * head, bool set_visited = false);

double get_rot_angle(Vector2d & v1, Vector2d & v2);

//make sure the min and max are in correct order
double intersect_two_box(Point2d lower_left, Point2d lower_right, Point2d upper_right,
                         Point2d upper_left, double xmin2, double xmax2, double ymin2,
                         double ymax2);

}