#include "RingClipper.h"
#include "intersection.h"
#include <map>
#include <fstream>
using namespace std;

namespace cusg
{
Point2d get_center(VNode * head)
{
    double x = 0, y = 0;
    int size = 0;
    VNode * node = head;
    do
    {
        x += node->x;
        y += node->y;
        size++;

        node = node->next;
    } while (node != head);

    return Point2d(x / size, y / size);
}

bool almost_same(double a[2], double b[2])
{
    return (fabs(a[0] - b[0]) <= SMALLCOORDINATE && fabs(a[1] - b[1]) <= SMALLCOORDINATE);
}
bool almost_equal(VNode * v1, VNode * v2)
{
    return (fabs(v1->x - v2->x) <= SMALLCOORDINATE && fabs(v1->y - v2->y) <= SMALLCOORDINATE);
}
//slope, along the line per unit distance from it, of the edge that leaves v on the kept side
double kept_slope(VNode * v, bool is_x)
{
    for (int side = 0; side < 2; side++)
    {
        VNode * n = side == 0 ? v->next : v->prev;
        while (n != v && almost_equal(n, v)) n = side == 0 ? n->next : n->prev;
        double along = is_x ? n->y - v->y : n->x - v->x;
        double away = fabs(is_x ? n->x - v->x : n->y - v->y);
        if (n->is_in_direct && away > SMALLCOORDINATE) return along / away;
    }
    return 0;
}
//true if a sorts after b along the line; crossings at the same place are sorted by where
//their rings run just inside the kept side
bool sorts_after(VNode * a, VNode * b, bool is_x)
{
    double pa = is_x ? a->y : a->x;
    double pb = is_x ? b->y : b->x;
    if (pa == pb || (almost_equal(a, b) && fabs(pa - pb) <= SMALLCOORDINATE))
        return kept_slope(a, is_x) > kept_slope(b, is_x);
    return pa > pb;
}
//sort the vertices based on X or Y ascendingly
void insertion_sort_v(vector<VNode *> & vertices, bool is_x)
{
    int vsize = vertices.size();
    if (vsize < 2) return;

    for (int i = 1; i < vsize; i++)
    {
        VNode * key = vertices[i];
        int j = i - 1;
        while (j >= 0 && sorts_after(vertices[j], key, is_x))
        {
            vertices[j + 1] = vertices[j];
            j--;
        }
        vertices[j + 1] = key;
    }
}

int locate_v_in_vector(vector<VNode *> & vertices, VNode * target, int start_index, int total_size)
{
    for (int i = start_index; i < total_size; i++)
    {
        if (vertices[i] == target) return i;
    }
    return -1;
}

void get_xy_bound(VNode * first, double & xmin, double & xmax, double & ymin, double & ymax)
{
    xmin = ymin = FLT_MAX;
    xmax = ymax = FLT_MIN;
    VNode * node = first;
    do
    {
        if (node->x < xmin) xmin = node->x;
        if (node->y < ymin) ymin = node->y;
        if (node->x > xmax) xmax = node->x;
        if (node->y > ymax) ymax = node->y;
        node = node->next;
    } while (node != first);
}

VNode * copy_ring(VNode * head, bool set_visited)
{
    if (head == nullptr) return nullptr;

    VNode * new_head = new VNode(head->x, head->y);
    if (set_visited) head->is_visited = true;
    new_head->next = new_head->prev = new_head;
    VNode * node = head->next;
    while (node != head)
    {
        VNode * node_copy = new VNode(node->x, node->y);
        new_head->prev->next = node_copy;
        node_copy->prev = new_head->prev;
        node_copy->next = new_head;
        new_head->prev = node_copy;

        if (set_visited) node->is_visited = true;

        node = node->next;
    }

    return new_head;
}
RingClipper * copy_ring_clipper(RingClipper * ring_clip, bool set_visited = false)
{
    RingClipper * clip_copy = new RingClipper();

    VNode * new_head = copy_ring(ring_clip->first, set_visited);
    clip_copy->first = new_head;
    clip_copy->rings.push_back(new_head);

    if (ring_clip->rings.size() > 1)
    {
        list<VNode *>::iterator vit = ring_clip->rings.begin();
        ++vit;
        for (; vit != ring_clip->rings.end(); ++vit)
        {
            VNode * ring_copy = copy_ring(*vit, set_visited);
            clip_copy->rings.push_back(ring_copy);
        }
    }
    return clip_copy;
}
void RingClipper::reset_status()
{
    //delete the new added
    for (auto added : added_vertices)
    {
        VNode * before = added->prev;
        VNode * after = added->next;
        before->next = after;
        after->prev = before;
        added->prev = added->next = nullptr;
        delete added;
    }
    added_vertices.clear();

    //loop over the influenced vertices
    for (auto vertex : changed_vertices)
    {
        vertex->reset_status();
    }
    changed_vertices.clear();
}

double get_rot_angle(Vector2d & v1, Vector2d & v2)
{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    //have bug here
    double angle = 0;
    Vector2d dir1 = v1.normalize();
    Vector2d dir2 = v2.normalize();
    double dot = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]);
    if (fabs(dot - 1.0) <= EPSILON)
        return 0;
    else if (fabs(dot + 1.0) <= EPSILON)
        return PI;
    else
    {
        double cross = 0;
        angle = acos(dot);
        cross = (dir1[0] * dir2[1]) - (dir2[0] * dir1[1]);
        if (cross < 0) angle = 2 * PI - angle;

        return angle;
    }
}
//test enclosing
bool clip_ring_encloses(VNode * head, const Point2d & p2d)
{
    if (head == nullptr || head->next == head) return false;
    //count the edges crossing a ray from the point; summing the angles they turn through
    //loses the small ones
    bool inside = false;
    VNode * node = head;
    do
    {
        double a[2] = {node->x, node->y};
        double b[2] = {node->next->x, node->next->y};
        if (crosses_ray(a, b, p2d.get())) inside = !inside;
        node = node->next;
    } while (node != head);

    return inside;
}

void get_feasible_pre_next(VNode * vertex, VNode *& prev_vertex, VNode *& next_vertex,
                           VNode *& prev_adjacent, VNode *& next_adjacent)
{
    prev_vertex = vertex->prev;
    next_vertex = vertex->next;
    while (almost_equal(prev_vertex, vertex))
    {
        prev_vertex = prev_vertex->prev;
    }
    while (almost_equal(next_vertex, vertex))
    {
        next_vertex = next_vertex->next;
    }
    prev_adjacent = prev_vertex;
    next_adjacent = next_vertex;
    //filtering
    while (prev_vertex->is_scan_point)
    {
        prev_adjacent = prev_vertex;
        prev_vertex = prev_vertex->prev;
    }
    while (next_vertex->is_scan_point)
    {
        next_adjacent = next_vertex;
        next_vertex = next_vertex->next;
    }
}
bool decide_next_order_for_out(VNode * vertex, bool is_x, bool is_bigger)
{
    assert(vertex->is_in_direct);
    assert(vertex->next->is_in_direct || vertex->prev->is_in_direct);

    if (!(vertex->prev->is_in_direct && vertex->next->is_in_direct))
        return vertex->next->is_in_direct;
    else
    {
        //1:pre; 2:next
        VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
              *next_adjacent = nullptr;
        get_feasible_pre_next(vertex, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
        assert(!(prev_adjacent->is_visited && next_adjacent->is_visited));
        //if one of them is not in indirection
        if (!(prev_vertex->is_in_direct && next_vertex->is_in_direct))
        {
            if (prev_vertex->is_in_direct && (!prev_adjacent->is_visited)) //previous
                return false;
            else if (next_vertex->is_in_direct && (!next_adjacent->is_visited)) //next
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (prev_vertex->is_in_direct && next_vertex->is_in_direct)
        {
            if (prev_adjacent->is_visited || prev_vertex->is_visited) return true;
            if (next_adjacent->is_visited || next_vertex->is_visited) return false;
        }
        if (!(prev_adjacent->is_scan_point && next_adjacent->is_scan_point))
        {
            Vector2d to_prev(prev_adjacent->x - vertex->x, prev_adjacent->y - vertex->y);
            Vector2d to_next(next_adjacent->x - vertex->x, next_adjacent->y - vertex->y);
            double angle = get_rot_angle(to_prev, to_next);
            if ((!is_x && is_bigger) || (is_x && !is_bigger))
                return (angle <= PI);
            else
                return (angle >= PI);
        }
        else
        {
            Vector2d to_prev(prev_adjacent->x - vertex->x, prev_adjacent->y - vertex->y);
            Vector2d to_next(next_vertex->x - vertex->x, next_vertex->y - vertex->y);
            double angle = get_rot_angle(to_prev, to_next);
            if ((!is_x && is_bigger) || (is_x && !is_bigger))
                return (angle <= PI);
            else
                return (angle >= PI);
        }
    }
}
bool decide_next_order_for_hole(VNode * vertex, bool is_x, bool is_bigger)
{
    assert(vertex->is_in_direct);
    assert(vertex->next->is_in_direct || vertex->prev->is_in_direct);

    if (!(vertex->prev->is_in_direct && vertex->next->is_in_direct))
        return vertex->next->is_in_direct;
    else
    {
        //1:pre; 2:next
        VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
              *next_adjacent = nullptr;
        get_feasible_pre_next(vertex, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
        assert(!(prev_adjacent->is_visited && next_adjacent->is_visited));
        if (!(prev_vertex->is_in_direct && next_vertex->is_in_direct))
        {
            if (prev_vertex->is_in_direct && (!prev_adjacent->is_visited))
                return false;
            else if (next_vertex->is_in_direct && (!next_adjacent->is_visited))
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (prev_vertex->is_in_direct && next_vertex->is_in_direct)
        {
            if (prev_adjacent->is_visited || prev_vertex->is_visited) return true;
            if (next_adjacent->is_visited || next_vertex->is_visited) return false;
        }
    }

    Vector2d to_prev(vertex->prev->x - vertex->x, vertex->prev->y - vertex->y);
    Vector2d to_next(vertex->next->x - vertex->x, vertex->next->y - vertex->y);

    double angle = get_rot_angle(to_prev, to_next);
    if (((!is_x) && is_bigger) || (is_x && (!is_bigger)))
        return (angle >= PI);
    else
        return (angle <= PI);
}

bool check_all_indirect(double line_coordinate, bool is_x, bool is_bigger, VNode * head)
{

    Point2d p2d = get_center(head);
    if (is_x)
    {
        return (is_bigger == (p2d[0] > line_coordinate));
    }
    else
    {
        return (is_bigger == (p2d[1] > line_coordinate));
    }
}
RingClipper::RingClipper()
{
    first = nullptr;
}
RingClipper::~RingClipper()
{
    destroy();
}
void RingClipper::destroy()
{
    for (auto & ring : rings)
    {
        destroy(ring);
    }
    rings.clear();
    first = nullptr;
}
void RingClipper::destroy(VNode * head)
{
    if (head == nullptr) return;
    VNode * node = head->next;
    while (node != head)
    {
        VNode * doomed = node;
        node = node->next;
        if (doomed != nullptr) delete doomed;
    }
    delete head;
    head = nullptr;
}

void RingClipper::add(VNode * v)
{
    //a repeated vertex is left out and freed
    if (first != nullptr && fabs(v->x - first->x) < SMALLNUMBER &&
        fabs(v->y - first->y) < SMALLNUMBER)
    {
        delete v;
        return;
    }
    if (first != nullptr && fabs(v->x - first->prev->x) < SMALLNUMBER &&
        fabs(v->y - first->prev->y) < SMALLNUMBER)
    {
        delete v;
        return;
    }

    if (first == nullptr)
    {
        first = v;
        first->next = v;
        first->prev = v;
    }
    else
    {
        VNode * next = first;
        VNode * prev = next->prev;
        next->prev = v;
        v->next = next;
        v->prev = prev;
        prev->next = v;
    }
}
void RingClipper::insert(VNode * before, VNode * after, VNode * v)
{

    before->next = v;
    v->prev = before;
    v->next = after;
    after->prev = v;
}
VNode * RingClipper::collect_along_way(vector<VNode *> & vertices, VNode * start, bool forward)
{
    assert(start->is_in_direct && start->is_scan_point);
    VNode * previous = start;
    vertices.push_back(previous);
    VNode * node = forward ? start->next : start->prev;

    while (almost_equal(previous, node) && node->is_in_direct)
    {
        vertices.push_back(node);
        previous = node;
        node = forward ? node->next : node->prev;
    }
    ///{
    while (node->is_in_direct)
    {
        vertices.push_back(node);
        if (node->is_scan_point) break;
        node = forward ? node->next : node->prev;
    }

    return vertices.back();
}

//true if v is on the line and both its neighbors are on the kept side, off the line
bool touches_kept_side(VNode * v)
{
    return v->is_online && v->next->is_in_direct && !v->next->is_online && v->prev->is_in_direct &&
           !v->prev->is_online;
}

//set intersection
int RingClipper::set_intersection(VNode * start_vertex, double line_coordinate, bool is_x,
                                  bool is_bigger, vector<VNode *> & line_vertices)
{
    double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
    get_xy_bound(start_vertex, xmin, xmax, ymin, ymax);
    double xwidth = fabs(xmax - xmin);
    double ywidth = fabs(ymax - ymin);
    if (xwidth < SMALLNUMBER) xwidth = ywidth;
    if (ywidth < SMALLNUMBER) ywidth = xwidth;

    double a[2] = {0}, b[2] = {0};
    if (is_x)
    {
        a[0] = line_coordinate;
        a[1] = ymin - ywidth / 2;
        b[0] = line_coordinate;
        b[1] = ymax + ywidth / 2;
    }
    else
    {
        a[0] = xmin - xwidth / 2;
        a[1] = line_coordinate;
        b[0] = xmax + xwidth / 2;
        b[1] = line_coordinate;
    }

    //move vertices that are on the line within the tolerance onto it, in place, so that an
    //edge along the line is intersected as one
    VNode * snap = start_vertex;
    do
    {
        if (snap->on_line(line_coordinate, is_x))
        {
            if (is_x)
                snap->x = line_coordinate;
            else
                snap->y = line_coordinate;
        }
        snap = snap->next;
    } while (snap != start_vertex);

    //try to intersect
    map<VNode *, bool> should_record;
    vector<VNode *> inter_cand;
    VNode * vertex = start_vertex;
    do
    {
        double c[2] = {0}, d[2] = {0}, p[2] = {FLT_MIN, FLT_MIN};
        vertex->set_direction(line_coordinate, is_x, is_bigger);
        //intersect
        c[0] = vertex->x;
        c[1] = vertex->y;
        d[0] = vertex->next->x;
        d[1] = vertex->next->y;
        //an edge with both ends on the line lies along it; seg_seg_int can miss this by rounding
        char code = 0;
        if ((is_x && c[0] == line_coordinate && d[0] == line_coordinate) ||
            (!is_x && c[1] == line_coordinate && d[1] == line_coordinate))
            code = 'e';
        else
            code = seg_seg_int(a, b, c, d, p);

        //set influenced vertices
        if (vertex->is_in_direct) changed_vertices.push_back(vertex);
        if (code == 'e')
        {
            vertex->is_online = vertex->next->is_online = true;
            vertex->is_scan_point = vertex->next->is_scan_point = true; //
            changed_vertices.push_back(vertex);
            changed_vertices.push_back(vertex->next);
        }
        else if (code == 'v')
        {
            //if(almost_equal(c, p))
            if (almost_same(c, p))
            {
                vertex->is_online = true;
                vertex->is_scan_point = true; //
                changed_vertices.push_back(vertex);
            }
            else
            {
                assert(almost_same(d, p));
                vertex->next->is_online = true;
                vertex->next->is_scan_point = true; //
                changed_vertices.push_back(vertex->next);
            }
        }
        else if (code == '1')
        {
            //if(almost_equal(c, p))
            if (almost_same(c, p))
            {
                vertex->is_online = true;
                vertex->is_scan_point = true; //
                changed_vertices.push_back(vertex);
            }
            //else if(almost_equal(d, p))
            else if (almost_same(d, p))
            {
                vertex->next->is_online = true;
                vertex->next->is_scan_point = true; //
                changed_vertices.push_back(vertex->next);
            }
            else
            {
                VNode * new_vertex = new VNode(p[0], p[1]);
                this->insert(vertex, vertex->next, new_vertex);
                //add to the added_vertices array
                added_vertices.push_back(new_vertex);

                new_vertex->is_in_direct = true;
                new_vertex->is_scan_point = true; //

                should_record.insert(make_pair(new_vertex, true));
                vertex = vertex->next;
            }
        }

        vertex = vertex->next;
    } while (vertex != start_vertex);

    //loop to decide online v
    vertex = start_vertex;
    do
    {
        if (vertex->is_online)
        {
            vertex->is_in_direct = true;
            inter_cand.push_back(vertex);
            should_record.insert(make_pair(vertex, true));
        }
        vertex = vertex->next;
    } while (vertex != start_vertex);

    vector<VNode *> useless_nodes;
    vector<VNode *> not_online;
    vector<VNode *> touching;
    for (auto candidate : inter_cand)
    {
        bool pure_in_dir1 = (candidate->next->is_in_direct && (!candidate->next->is_online));
        bool pure_in_dir2 = (candidate->prev->is_in_direct && (!candidate->prev->is_online));
        if (candidate->is_online && (!(pure_in_dir1 && pure_in_dir2)))
        {
            not_online.push_back(candidate);
        }

        //a candidate that touches the line from the kept side does not cross it, unless the
        //ring also crosses the line at the same place
        if (touches_kept_side(candidate))
        {
            bool alone = true;
            for (auto & mit : should_record)
            {
                if (mit.first != candidate && almost_equal(mit.first, candidate) &&
                    !touches_kept_side(mit.first))
                    alone = false;
            }
            if (alone)
            {
                touching.push_back(candidate);
                should_record[candidate] = false;
                continue;
            }
        }

        if (pure_in_dir1 || pure_in_dir2)
            should_record[candidate] = true;
        else
            should_record[candidate] = false;

        if (!(pure_in_dir1 || pure_in_dir2))
        {
            useless_nodes.push_back(candidate);
        }
    }
    //not online
    for (auto candidate : not_online)
    {
        candidate->is_online = false;
    }
    //useless nodes
    for (auto candidate : useless_nodes)
    {
        candidate->is_online = false;
        candidate->is_in_direct = false; // candidate->is_scan_point = false;
    }
    //touching nodes are kept as ordinary vertices
    for (auto candidate : touching)
    {
        candidate->is_online = false;
        candidate->is_scan_point = false;
    }

    //decide the final intersection vertices, in ring order rather than by address so that
    //vertices at the same place on the scanning line are always sorted the same way
    int on_num = 0;
    vertex = start_vertex;
    do
    {
        map<VNode *, bool>::iterator mit = should_record.find(vertex);
        if (mit != should_record.end() && mit->second)
        {
            line_vertices.push_back(vertex);
            on_num++;
        }
        vertex = vertex->next;
    } while (vertex != start_vertex);

    return on_num;
}

void RingClipper::collect_second_part(vector<VNode *> & line_vertices, int start_index,
                                      int end_index, VNode * start_node, VNode * end_node,
                                      bool forward, int vertex_count,
                                      list<RingClipper *> & clip_results, bool is_x, bool is_bigger,
                                      vector<VNode *> & second_half)
{

    //collect these vertices
    vector<VNode *> along_vertices;
    VNode * along_end = collect_along_way(along_vertices, line_vertices[start_index], forward);

    for (auto vertex : along_vertices)
    {
        if (!vertex->is_visited)
        {
            second_half.push_back(vertex);
            vertex->is_visited = true;
        }
    }

    int along_end_index = locate_v_in_vector(line_vertices, along_end, start_index, vertex_count);
    assert(along_end_index != -1);

    /**************************************************************/
    //deal with the first part
    int first_start = start_index + 1;
    bool first_forward =
        false; //decide_next_order_for_out(line_vertices[first_start], is_x, is_bigger);
    VNode * start_neighbor =
        forward ? line_vertices[start_index]->prev : line_vertices[start_index]->next;
    if (line_vertices[start_index]->is_online && start_neighbor->is_in_direct &&
        (!start_neighbor->is_visited))
    {
        first_start = start_index;
        first_forward = !forward;
    }
    if (first_start != start_index && first_start < along_end_index)
        first_forward = decide_next_order_for_out(line_vertices[first_start], is_x, is_bigger);

    int first_end = along_end_index - 1;
    VNode * after_end_first =
        forward ? line_vertices[along_end_index]->next : line_vertices[along_end_index]->prev;
    if (line_vertices[along_end_index]->is_online && after_end_first->is_in_direct &&
        (!after_end_first->is_visited))
    {
        bool end_forward = decide_next_order_for_hole(
            line_vertices[along_end_index], is_x,
            is_bigger); // !(decide_next_order_for_out(line_vertices[along_end_index], is_x, is_bigger));
        VNode * hole_next = end_forward ? line_vertices[along_end_index]->next
                                        : line_vertices[along_end_index]->prev;
        if (hole_next->is_visited)
        {
            first_end = along_end_index;
        }
    }

    //deal with the first part
    if (first_start < first_end)
    {
        collect_polygons_in_interval(is_x, is_bigger, line_vertices, clip_results, first_start,
                                     first_end, first_forward);
    }
    //end of first part
    /**********************************************************************************************/

    /*********************************************************************************************/
    //deal with the second part
    int second_start = along_end_index + 1;
    bool second_forward =
        false; //decide_next_order_for_out(line_vertices[second_start], is_x, is_bigger);
    VNode * after_end_second =
        forward ? line_vertices[along_end_index]->next : line_vertices[along_end_index]->prev;
    if (line_vertices[along_end_index]->is_online && after_end_second->is_in_direct &&
        (!after_end_second->is_visited))
    {
        bool end_forward = decide_next_order_for_hole(
            line_vertices[along_end_index], is_x,
            is_bigger); //!(decide_next_order_for_out(line_vertices[along_end_index], is_x, is_bigger));
        VNode * hole_next = end_forward ? line_vertices[along_end_index]->next
                                        : line_vertices[along_end_index]->prev;
        if (!hole_next->is_visited)
        {
            second_start = along_end_index;
            second_forward = forward;
        }
    }

    if (second_start != along_end_index && second_start < end_index)
    {
        second_forward = decide_next_order_for_out(line_vertices[second_start], is_x, is_bigger);
    }

    if (second_start < end_index)
    {
        collect_second_part(line_vertices, second_start, end_index, line_vertices[second_start],
                            end_node, second_forward, vertex_count, clip_results, is_x, is_bigger,
                            second_half);
    }
    //end of second part
    /*********************************************************************************************/
}

void RingClipper::construct_clip_rings(vector<VNode *> & first_half, vector<VNode *> & second_half,
                                       list<RingClipper *> & clip_results)
{
    if (first_half.size() + second_half.size() < 3)
        return;
    else
    {
        RingClipper * ring_clip = new RingClipper();
        for (auto & vit : first_half)
        {
            ring_clip->add(copy(vit));
        }
        for (vector<VNode *>::reverse_iterator rvit = second_half.rbegin();
             rvit != second_half.rend(); ++rvit)
        {
            ring_clip->add(copy(*rvit));
        }
        ring_clip->rings.push_back(ring_clip->first);

        Ring ring;
        convert_clip_to_ring(ring_clip->first, ring);
        if (fabs(ring.get_area()) < SMALLAREA) //SMALLNUMBER
        {
            ring.destroy();
            delete ring_clip;
            return;
        }
        else
        {
            ring.destroy();
            clip_results.push_back(ring_clip);
        }
    }
}
void RingClipper::collect_ring(vector<VNode *> & line_vertices, int start_index, int end_index,
                               bool is_next_order, list<RingClipper *> & clip_results, bool is_x,
                               bool is_bigger)
{
    int vsize = line_vertices.size();
    vector<VNode *> first_half;
    vector<VNode *> second_half;

    //collect first half: out boundary
    VNode * start_node = line_vertices[start_index];
    VNode * end_node = line_vertices[end_index];
    VNode * node = start_node;
    first_half.push_back(start_node);
    start_node->is_visited = true;
    do
    {
        if (is_next_order)
            node = node->next;
        else
            node = node->prev;

        node->is_visited = true;
        first_half.push_back(node);
    } while (node != end_node);

    //filter the starting vertex
    int start_id = start_index + 1;
    int end_id = end_index - 1;
    bool forward = !is_next_order;
    VNode * start_neighbor =
        !is_next_order ? line_vertices[start_index]->next : line_vertices[start_index]->prev;
    if (line_vertices[start_index]->is_online && (!start_neighbor->is_visited))
    {
        assert(start_neighbor->is_in_direct);
        forward = !is_next_order;
        start_id = start_index;
    }
    else
    {
        start_id = start_index + 1;
        if (start_id < end_id)
            forward = decide_next_order_for_out(line_vertices[start_id], is_x, is_bigger);
    }

    //if(line_vertices[end_index]->is_online

    if (start_id < end_id)
    {
        collect_second_part(line_vertices, start_id, end_id, line_vertices[start_id],
                            line_vertices[end_id], forward, vsize, clip_results, is_x, is_bigger,
                            second_half);
    }

    //combine two lists to get a new polygon
    construct_clip_rings(first_half, second_half, clip_results);
}

void RingClipper::collect_polygons_in_interval(bool is_x, bool is_bigger,
                                               vector<VNode *> & line_vertices,
                                               list<RingClipper *> & clips, int from, int to,
                                               bool start_forward)
{
    int size = line_vertices.size();
    bool forward = start_forward;

    int start_id = from;
    while (start_id < to)
    {
        VNode * start_node = line_vertices[start_id];
        assert(start_node->next->is_in_direct || start_node->prev->is_in_direct);

        vector<VNode *> along_vertices;
        VNode * end_node = collect_along_way(along_vertices, start_node, forward);
        int end_id = locate_v_in_vector(line_vertices, end_node, start_id, size);

        assert(end_id != -1);

        collect_ring(line_vertices, start_id, end_id, forward, clips, is_x,
                     is_bigger); //note here or not, may have bugs here

        if (end_id >= to)
        {
            assert(end_id == to);
            break;
        }

        VNode * end_onward = forward ? line_vertices[end_id]->next : line_vertices[end_id]->prev;
        if (line_vertices[end_id]->is_online && (!end_onward->is_visited))
        {
            start_id = end_id;
        }
        else
        {
            start_id = end_id + 1;
            forward = decide_next_order_for_out(line_vertices[start_id], is_x, is_bigger);
        }
    }
}

void RingClipper::collect_polygons(double line_coordinate, bool is_x, bool is_bigger,
                                   vector<int> & line_vertex_counts,
                                   vector<VNode *> & line_vertices, list<RingClipper *> & clips)
{
    int size = line_vertices.size();
    if (size <= 1)
    {
        if (check_all_indirect(line_coordinate, is_x, is_bigger, first))
        {
            clips.push_back(copy_ring_clipper(this, true));
        }
        //set visited
        //for(list<VNode*>::iterator )
        return;
    }

    //sort the intersection vertices
    insertion_sort_v(line_vertices, is_x);

    bool is_next_order =
        decide_next_order_for_out(line_vertices[0], is_x, is_bigger); //need modify here

    collect_polygons_in_interval(is_x, is_bigger, line_vertices, clips, 0, size - 1, is_next_order);
}

void RingClipper::clip(double line_coordinate, bool is_x, bool is_bigger,
                       list<RingClipper *> & clips)
{
    vector<int> line_vertex_counts;
    vector<VNode *> line_vertices;
    //set intersections for holes
    for (auto ring : rings)
    {
        int line_vertex_count =
            set_intersection(ring, line_coordinate, is_x, is_bigger, line_vertices);
        line_vertex_counts.push_back(line_vertex_count);
    }

    //collect the out boundary
    collect_polygons(line_coordinate, is_x, is_bigger, line_vertex_counts, line_vertices, clips);

    /**********************************************/
    if (rings.size() == 0 || line_vertex_counts.size() == 0) return;

    list<VNode *>::iterator pit = rings.begin();
    vector<int>::iterator iit = line_vertex_counts.begin();
    //++pit; ++iit; //discard the out boundary
    for (; pit != rings.end(); ++pit, ++iit)
    {

        if (*iit == 0 && ((*pit)->is_visited == false) && (*pit)->is_in_direct)
        {
            if ((*pit)->is_scan_point)
            {
                VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
                      *next_adjacent = nullptr;
                get_feasible_pre_next(*pit, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
                if (!prev_adjacent->is_in_direct || !next_adjacent->is_in_direct) continue;
            }

            bool have_added = false;
            VNode * new_vertex = *pit;
            Ring nrings(Ring::PIN);
            convert_clip_to_ring(new_vertex, nrings);
            const Point2d & center = nrings.get_center();
            for (auto & clip_part : clips)
            {
                if (clip_ring_encloses(clip_part->first, center))
                {
                    clip_part->rings.push_back(copy_ring(new_vertex));
                    have_added = true;
                    break;
                }
            }
            if (!have_added)
            {
                RingClipper * ring_clip = convert_ring_to_clip(nrings);
                ring_clip->rings.push_back(ring_clip->first);
                clips.push_back(ring_clip);
            }
            nrings.destroy(); //clear memory
        }
    }
}
void RingClipper::clip(double xmin, double xmax, double ymin, double ymax,
                       list<RingClipper *> & clips)
{
    //pre-processing
    list<RingClipper *> inside_xmin, inside_xmax, inside_ymin;
    this->clip(xmin, true, true, inside_xmin);

    for (auto ring_clip : inside_xmin)
    {
        ring_clip->clip(xmax, true, false, inside_xmax);

        ring_clip->destroy();
        delete ring_clip;
    }

    for (auto ring_clip : inside_xmax)
    {
        ring_clip->clip(ymin, false, true, inside_ymin);

        ring_clip->destroy();
        delete ring_clip;
    }

    for (auto ring_clip : inside_ymin)
    {
        ring_clip->clip(ymax, false, false, clips);

        ring_clip->destroy();
        delete ring_clip;
    }
}

VNode * RingClipper::copy(const VNode * v)
{
    VNode * new_vertex = new VNode(v->x, v->y);
    return new_vertex;
}

VNode * convert_ring_to_vertex_nodes(const Ring & ring)
{
    LinkedVertex * head = ring.get_head();
    VNode * new_head = new VNode(head->get_pos()[0], head->get_pos()[1]);
    new_head->next = new_head->prev = new_head;
    LinkedVertex * vertex = head->get_next();
    while (vertex != head)
    {
        VNode * node_copy = new VNode(vertex->get_pos()[0], vertex->get_pos()[1]);
        //leave out a repeated vertex, such as the closing vertex of a shapefile ring
        if (almost_equal(node_copy, new_head->prev) || almost_equal(node_copy, new_head))
        {
            delete node_copy;
            vertex = vertex->get_next();
            continue;
        }
        new_head->prev->next = node_copy;
        node_copy->prev = new_head->prev;
        node_copy->next = new_head;
        new_head->prev = node_copy;

        vertex = vertex->get_next();
    }

    return new_head;
}

//convert the ring to RingClipper
RingClipper * convert_ring_to_clip(const Ring & ring)
{
    RingClipper * ring_clip = new RingClipper();
    LinkedVertex * head = ring.get_head();
    LinkedVertex * vertex = head;
    do
    {
        VNode * node = new VNode(vertex->get_pos()[0], vertex->get_pos()[1]);
        ring_clip->add(node);

        vertex = vertex->get_next();
    } while (vertex != head);
    return ring_clip;
}

//convert
void convert_clip_to_ring(VNode * first_vertex, Ring & ring)
{
    ring.begin_poly();
    VNode * vertex = first_vertex;
    do
    {
        ring.add_vertex(vertex->x, vertex->y, true);
        vertex = vertex->next;
    } while (vertex != first_vertex);

    ring.end_poly(false);
}

//convert
RingClipper * convert_polygon_to_clip(const GisPolygon & polygon)
{
    RingClipper * ring_clip = new RingClipper();
    for (const auto & part : polygon)
    {
        if (part.get_type() == Ring::POUT)
        {
            VNode * part_vertices = convert_ring_to_vertex_nodes(part);
            ring_clip->first = part_vertices;
            ring_clip->rings.push_front(part_vertices);
        }
        else
        {
            VNode * part_vertices = convert_ring_to_vertex_nodes(part);
            //a hole with fewer than three distinct vertices has no area
            if (part_vertices->next == part_vertices || part_vertices->next->next == part_vertices)
                ring_clip->destroy(part_vertices);
            else
                ring_clip->rings.push_back(part_vertices);
        }
    }
    return ring_clip;
}
//convert
void convert_clip_to_polygon(RingClipper * ring_clip, GisPolygon * polygon)
{
    if (ring_clip == nullptr || ring_clip->first == nullptr) return;
    //loop over
    Ring outer_ring(Ring::POUT);
    convert_clip_to_ring(ring_clip->first, outer_ring);
    polygon->push_back(outer_ring);

    if (ring_clip->rings.size() == 1) return;
    list<VNode *>::iterator vit = ring_clip->rings.begin();
    ++vit;
    for (; vit != ring_clip->rings.end(); ++vit)
    {
        Ring hole(Ring::PIN);
        convert_clip_to_ring(*vit, hole);
        polygon->push_back(hole);
    }
}

RingClipper * convert_simple_polygon_to_clip(const Polygon & polygon)
{
    RingClipper * ring_clip = new RingClipper();
    for (const auto & part : polygon)
    {
        if (part.get_type() == Ring::POUT)
        {
            VNode * part_vertices = convert_ring_to_vertex_nodes(part);
            ring_clip->first = part_vertices;
            ring_clip->rings.push_front(part_vertices);
        }
        else
        {
            VNode * part_vertices = convert_ring_to_vertex_nodes(part);
            //a hole with fewer than three distinct vertices has no area
            if (part_vertices->next == part_vertices || part_vertices->next->next == part_vertices)
                ring_clip->destroy(part_vertices);
            else
                ring_clip->rings.push_back(part_vertices);
        }
    }
    return ring_clip;
}

VNode * convert_box_to_clip(Point2d lower_left, Point2d lower_right, Point2d upper_right,
                            Point2d upper_left)
{
    VNode * corners[4];
    corners[0] = new VNode(lower_left[0], lower_left[1]);
    corners[1] = new VNode(lower_right[0], lower_right[1]);
    corners[2] = new VNode(upper_right[0], upper_right[1]);
    corners[3] = new VNode(upper_left[0], upper_left[1]);

    VNode * head = corners[0];
    head->next = head->prev = corners[0];

    for (int i = 1; i <= 3; i++)
    {
        head->prev->next = corners[i];
        corners[i]->prev = head->prev;
        corners[i]->next = head;
        head->prev = corners[i];
    }

    return head;
}

double get_out_bd_area(RingClipper * ring_clip)
{
    Ring outer_ring(Ring::POUT);
    convert_clip_to_ring(ring_clip->first, outer_ring);

    double area = outer_ring.get_area();
    outer_ring.destroy();
    return area;
}

//make sure the min and max are in correct order
double intersect_two_box(Point2d lower_left, Point2d lower_right, Point2d upper_right,
                         Point2d upper_left, double xmin2, double xmax2, double ymin2, double ymax2)
{
    double area(0.0);

    //build a BoxClip
    RingClipper ring_clip;
    VNode * head = convert_box_to_clip(lower_left, lower_right, upper_right, upper_left);
    ring_clip.first = head;
    ring_clip.rings.push_back(head);

    list<RingClipper *> clips;
    ring_clip.clip(xmin2, xmax2, ymin2, ymax2, clips);

    for (auto & clip_part : clips)
    {

        area += get_out_bd_area(clip_part);

        clip_part->destroy();
        delete clip_part;
    }

    ring_clip.destroy();
    return area;
}

}
