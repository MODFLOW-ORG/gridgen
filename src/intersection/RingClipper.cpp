#include "RingClipper.h"
#include "intersection.h"
#include <map>
#include <fstream>
using namespace std;

namespace cusg
{
Point2d getCenter(VNode * head)
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

bool AlmostSame(double a[2], double b[2])
{
    return (fabs(a[0] - b[0]) <= SMALLCOORDINATE && fabs(a[1] - b[1]) <= SMALLCOORDINATE);
}
bool AlmostEqual(VNode * v1, VNode * v2)
{
    return (fabs(v1->x - v2->x) <= SMALLCOORDINATE && fabs(v1->y - v2->y) <= SMALLCOORDINATE);
}
//slope, along the line per unit distance from it, of the edge that leaves v on the kept side
double keptSlope(VNode * v, bool isX)
{
    for (int side = 0; side < 2; side++)
    {
        VNode * n = side == 0 ? v->next : v->prev;
        while (n != v && AlmostEqual(n, v)) n = side == 0 ? n->next : n->prev;
        double along = isX ? n->y - v->y : n->x - v->x;
        double away = fabs(isX ? n->x - v->x : n->y - v->y);
        if (n->isInDirect && away > SMALLCOORDINATE) return along / away;
    }
    return 0;
}
//true if a sorts after b along the line; crossings at the same place are sorted by where
//their rings run just inside the kept side
bool sortsAfter(VNode * a, VNode * b, bool isX)
{
    double pa = isX ? a->y : a->x;
    double pb = isX ? b->y : b->x;
    if (pa == pb || (AlmostEqual(a, b) && fabs(pa - pb) <= SMALLCOORDINATE))
        return keptSlope(a, isX) > keptSlope(b, isX);
    return pa > pb;
}
//sort the vertices based on X or Y ascendingly
void InsertionSortV(vector<VNode *> & vertices, bool isX)
{
    int vsize = vertices.size();
    if (vsize < 2) return;

    for (int i = 1; i < vsize; i++)
    {
        VNode * key = vertices[i];
        int j = i - 1;
        while (j >= 0 && sortsAfter(vertices[j], key, isX))
        {
            vertices[j + 1] = vertices[j];
            j--;
        }
        vertices[j + 1] = key;
    }
}

int locateVInVector(vector<VNode *> & vertices, VNode * target, int start_index, int totalSize)
{
    for (int i = start_index; i < totalSize; i++)
    {
        if (vertices[i] == target) return i;
    }
    return -1;
}

void getXYBound(VNode * first, double & xmin, double & xmax, double & ymin, double & ymax)
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

VNode * copy_ring(VNode * head, bool bSetVisited)
{
    if (head == nullptr) return nullptr;

    VNode * new_head = new VNode(head->x, head->y);
    if (bSetVisited) head->isVisited = true;
    new_head->next = new_head->prev = new_head;
    VNode * node = head->next;
    while (node != head)
    {
        VNode * node_copy = new VNode(node->x, node->y);
        new_head->prev->next = node_copy;
        node_copy->prev = new_head->prev;
        node_copy->next = new_head;
        new_head->prev = node_copy;

        if (bSetVisited) node->isVisited = true;

        node = node->next;
    }

    return new_head;
}
RingClipper * copy_ring_clipper(RingClipper * ring_clip, bool bSetVisited = false)
{
    RingClipper * clip_copy = new RingClipper();

    VNode * new_head = copy_ring(ring_clip->first, bSetVisited);
    clip_copy->first = new_head;
    clip_copy->rings.push_back(new_head);

    if (ring_clip->rings.size() > 1)
    {
        list<VNode *>::iterator vit = ring_clip->rings.begin();
        ++vit;
        for (; vit != ring_clip->rings.end(); ++vit)
        {
            VNode * ring_copy = copy_ring(*vit, bSetVisited);
            clip_copy->rings.push_back(ring_copy);
        }
    }
    return clip_copy;
}
void RingClipper::resetStatus()
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
        vertex->resetStatus();
    }
    changed_vertices.clear();
}

double getRotAngle(Vector2d & v1, Vector2d & v2)
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
        if (CrossesRay(a, b, p2d.get())) inside = !inside;
        node = node->next;
    } while (node != head);

    return inside;
}

void getFeasiblePreNext(VNode * vertex, VNode *& prev_vertex, VNode *& next_vertex,
                        VNode *& prev_adjacent, VNode *& next_adjacent)
{
    prev_vertex = vertex->prev;
    next_vertex = vertex->next;
    while (AlmostEqual(prev_vertex, vertex))
    {
        prev_vertex = prev_vertex->prev;
    }
    while (AlmostEqual(next_vertex, vertex))
    {
        next_vertex = next_vertex->next;
    }
    prev_adjacent = prev_vertex;
    next_adjacent = next_vertex;
    //filtering
    while (prev_vertex->isScanPoint)
    {
        prev_adjacent = prev_vertex;
        prev_vertex = prev_vertex->prev;
    }
    while (next_vertex->isScanPoint)
    {
        next_adjacent = next_vertex;
        next_vertex = next_vertex->next;
    }
}
bool decideNextOrderForOut(VNode * vertex, bool isX, bool isBigger)
{
    assert(vertex->isInDirect);
    assert(vertex->next->isInDirect || vertex->prev->isInDirect);

    if (!(vertex->prev->isInDirect && vertex->next->isInDirect))
        return vertex->next->isInDirect;
    else
    {
        //1:pre; 2:next
        VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
              *next_adjacent = nullptr;
        getFeasiblePreNext(vertex, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
        assert(!(prev_adjacent->isVisited && next_adjacent->isVisited));
        //if one of them is not in indirection
        if (!(prev_vertex->isInDirect && next_vertex->isInDirect))
        {
            if (prev_vertex->isInDirect && (!prev_adjacent->isVisited)) //previous
                return false;
            else if (next_vertex->isInDirect && (!next_adjacent->isVisited)) //next
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (prev_vertex->isInDirect && next_vertex->isInDirect)
        {
            if (prev_adjacent->isVisited || prev_vertex->isVisited) return true;
            if (next_adjacent->isVisited || next_vertex->isVisited) return false;
        }
        if (!(prev_adjacent->isScanPoint && next_adjacent->isScanPoint))
        {
            Vector2d to_prev(prev_adjacent->x - vertex->x, prev_adjacent->y - vertex->y);
            Vector2d to_next(next_adjacent->x - vertex->x, next_adjacent->y - vertex->y);
            double angle = getRotAngle(to_prev, to_next);
            if ((!isX && isBigger) || (isX && !isBigger))
                return (angle <= PI);
            else
                return (angle >= PI);
        }
        else
        {
            Vector2d to_prev(prev_adjacent->x - vertex->x, prev_adjacent->y - vertex->y);
            Vector2d to_next(next_vertex->x - vertex->x, next_vertex->y - vertex->y);
            double angle = getRotAngle(to_prev, to_next);
            if ((!isX && isBigger) || (isX && !isBigger))
                return (angle <= PI);
            else
                return (angle >= PI);
        }
    }
}
bool decideNextOrderForHole(VNode * vertex, bool isX, bool isBigger)
{
    assert(vertex->isInDirect);
    assert(vertex->next->isInDirect || vertex->prev->isInDirect);

    if (!(vertex->prev->isInDirect && vertex->next->isInDirect))
        return vertex->next->isInDirect;
    else
    {
        //1:pre; 2:next
        VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
              *next_adjacent = nullptr;
        getFeasiblePreNext(vertex, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
        assert(!(prev_adjacent->isVisited && next_adjacent->isVisited));
        if (!(prev_vertex->isInDirect && next_vertex->isInDirect))
        {
            if (prev_vertex->isInDirect && (!prev_adjacent->isVisited))
                return false;
            else if (next_vertex->isInDirect && (!next_adjacent->isVisited))
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (prev_vertex->isInDirect && next_vertex->isInDirect)
        {
            if (prev_adjacent->isVisited || prev_vertex->isVisited) return true;
            if (next_adjacent->isVisited || next_vertex->isVisited) return false;
        }
    }

    Vector2d to_prev(vertex->prev->x - vertex->x, vertex->prev->y - vertex->y);
    Vector2d to_next(vertex->next->x - vertex->x, vertex->next->y - vertex->y);

    double angle = getRotAngle(to_prev, to_next);
    if (((!isX) && isBigger) || (isX && (!isBigger)))
        return (angle >= PI);
    else
        return (angle <= PI);
}

bool checkAllIndirect(double line_coordinate, bool isX, bool isBigger, VNode * head)
{

    Point2d p2d = getCenter(head);
    if (isX)
    {
        return (isBigger == (p2d[0] > line_coordinate));
    }
    else
    {
        return (isBigger == (p2d[1] > line_coordinate));
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
VNode * RingClipper::collectAlongWay(vector<VNode *> & vertices, VNode * start, bool forward)
{
    assert(start->isInDirect && start->isScanPoint);
    VNode * previous = start;
    vertices.push_back(previous);
    VNode * node = forward ? start->next : start->prev;

    while (AlmostEqual(previous, node) && node->isInDirect)
    {
        vertices.push_back(node);
        previous = node;
        node = forward ? node->next : node->prev;
    }
    ///{
    while (node->isInDirect)
    {
        vertices.push_back(node);
        if (node->isScanPoint) break;
        node = forward ? node->next : node->prev;
    }

    return vertices.back();
}

//true if v is on the line and both its neighbors are on the kept side, off the line
bool touchesKeptSide(VNode * v)
{
    return v->isOnline && v->next->isInDirect && !v->next->isOnline && v->prev->isInDirect &&
           !v->prev->isOnline;
}

//set intersection
int RingClipper::setIntersection(VNode * start_vertex, double line_coordinate, bool isX,
                                 bool isBigger, vector<VNode *> & line_vertices)
{
    double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
    getXYBound(start_vertex, xmin, xmax, ymin, ymax);
    double xwidth = fabs(xmax - xmin);
    double ywidth = fabs(ymax - ymin);
    if (xwidth < SMALLNUMBER) xwidth = ywidth;
    if (ywidth < SMALLNUMBER) ywidth = xwidth;

    double a[2] = {0}, b[2] = {0};
    if (isX)
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
        if (snap->Online(line_coordinate, isX))
        {
            if (isX)
                snap->x = line_coordinate;
            else
                snap->y = line_coordinate;
        }
        snap = snap->next;
    } while (snap != start_vertex);

    //try to intersect
    map<VNode *, bool> shouldRecord;
    vector<VNode *> interCand;
    VNode * vertex = start_vertex;
    do
    {
        double c[2] = {0}, d[2] = {0}, p[2] = {FLT_MIN, FLT_MIN};
        vertex->setDirection(line_coordinate, isX, isBigger);
        //intersect
        c[0] = vertex->x;
        c[1] = vertex->y;
        d[0] = vertex->next->x;
        d[1] = vertex->next->y;
        //an edge with both ends on the line lies along it; SegSegInt can miss this by rounding
        char code = 0;
        if ((isX && c[0] == line_coordinate && d[0] == line_coordinate) ||
            (!isX && c[1] == line_coordinate && d[1] == line_coordinate))
            code = 'e';
        else
            code = SegSegInt(a, b, c, d, p);

        //set influenced vertices
        if (vertex->isInDirect) changed_vertices.push_back(vertex);
        if (code == 'e')
        {
            vertex->isOnline = vertex->next->isOnline = true;
            vertex->isScanPoint = vertex->next->isScanPoint = true; //
            changed_vertices.push_back(vertex);
            changed_vertices.push_back(vertex->next);
        }
        else if (code == 'v')
        {
            //if(AlmostEqual(c, p))
            if (AlmostSame(c, p))
            {
                vertex->isOnline = true;
                vertex->isScanPoint = true; //
                changed_vertices.push_back(vertex);
            }
            else
            {
                assert(AlmostSame(d, p));
                vertex->next->isOnline = true;
                vertex->next->isScanPoint = true; //
                changed_vertices.push_back(vertex->next);
            }
        }
        else if (code == '1')
        {
            //if(AlmostEqual(c, p))
            if (AlmostSame(c, p))
            {
                vertex->isOnline = true;
                vertex->isScanPoint = true; //
                changed_vertices.push_back(vertex);
            }
            //else if(AlmostEqual(d, p))
            else if (AlmostSame(d, p))
            {
                vertex->next->isOnline = true;
                vertex->next->isScanPoint = true; //
                changed_vertices.push_back(vertex->next);
            }
            else
            {
                VNode * new_vertex = new VNode(p[0], p[1]);
                this->insert(vertex, vertex->next, new_vertex);
                //add to the added_vertices array
                added_vertices.push_back(new_vertex);

                new_vertex->isInDirect = true;
                new_vertex->isScanPoint = true; //

                shouldRecord.insert(make_pair(new_vertex, true));
                vertex = vertex->next;
            }
        }

        vertex = vertex->next;
    } while (vertex != start_vertex);

    //loop to decide online v
    vertex = start_vertex;
    do
    {
        if (vertex->isOnline)
        {
            vertex->isInDirect = true;
            interCand.push_back(vertex);
            shouldRecord.insert(make_pair(vertex, true));
        }
        vertex = vertex->next;
    } while (vertex != start_vertex);

    vector<VNode *> uselessNodes;
    vector<VNode *> notOnline;
    vector<VNode *> touching;
    for (auto candidate : interCand)
    {
        bool pureInDir1 = (candidate->next->isInDirect && (!candidate->next->isOnline));
        bool pureInDir2 = (candidate->prev->isInDirect && (!candidate->prev->isOnline));
        if (candidate->isOnline && (!(pureInDir1 && pureInDir2)))
        {
            notOnline.push_back(candidate);
        }

        //a candidate that touches the line from the kept side does not cross it, unless the
        //ring also crosses the line at the same place
        if (touchesKeptSide(candidate))
        {
            bool alone = true;
            for (auto & mit : shouldRecord)
            {
                if (mit.first != candidate && AlmostEqual(mit.first, candidate) &&
                    !touchesKeptSide(mit.first))
                    alone = false;
            }
            if (alone)
            {
                touching.push_back(candidate);
                shouldRecord[candidate] = false;
                continue;
            }
        }

        if (pureInDir1 || pureInDir2)
            shouldRecord[candidate] = true;
        else
            shouldRecord[candidate] = false;

        if (!(pureInDir1 || pureInDir2))
        {
            uselessNodes.push_back(candidate);
        }
    }
    //not online
    for (auto candidate : notOnline)
    {
        candidate->isOnline = false;
    }
    //useless nodes
    for (auto candidate : uselessNodes)
    {
        candidate->isOnline = false;
        candidate->isInDirect = false; // candidate->isScanPoint = false;
    }
    //touching nodes are kept as ordinary vertices
    for (auto candidate : touching)
    {
        candidate->isOnline = false;
        candidate->isScanPoint = false;
    }

    //decide the final intersection vertices, in ring order rather than by address so that
    //vertices at the same place on the scanning line are always sorted the same way
    int onNum = 0;
    vertex = start_vertex;
    do
    {
        map<VNode *, bool>::iterator mit = shouldRecord.find(vertex);
        if (mit != shouldRecord.end() && mit->second)
        {
            line_vertices.push_back(vertex);
            onNum++;
        }
        vertex = vertex->next;
    } while (vertex != start_vertex);

    return onNum;
}

void RingClipper::collectSecondPart(vector<VNode *> & line_vertices, int start_index, int end_index,
                                    VNode * start_node, VNode * end_node, bool forward,
                                    int vertex_count, list<RingClipper *> & clip_results, bool isX,
                                    bool isBigger, vector<VNode *> & secondHalf)
{

    //collect these vertices
    vector<VNode *> along_vertices;
    VNode * along_end = collectAlongWay(along_vertices, line_vertices[start_index], forward);

    for (auto vertex : along_vertices)
    {
        if (!vertex->isVisited)
        {
            secondHalf.push_back(vertex);
            vertex->isVisited = true;
        }
    }

    int along_end_index = locateVInVector(line_vertices, along_end, start_index, vertex_count);
    assert(along_end_index != -1);

    /**************************************************************/
    //deal with the first part
    int first_start = start_index + 1;
    bool first_forward = false; //decideNextOrderForOut(line_vertices[first_start], isX, isBigger);
    VNode * start_neighbor =
        forward ? line_vertices[start_index]->prev : line_vertices[start_index]->next;
    if (line_vertices[start_index]->isOnline && start_neighbor->isInDirect &&
        (!start_neighbor->isVisited))
    {
        first_start = start_index;
        first_forward = !forward;
    }
    if (first_start != start_index && first_start < along_end_index)
        first_forward = decideNextOrderForOut(line_vertices[first_start], isX, isBigger);

    int first_end = along_end_index - 1;
    VNode * after_end_first =
        forward ? line_vertices[along_end_index]->next : line_vertices[along_end_index]->prev;
    if (line_vertices[along_end_index]->isOnline && after_end_first->isInDirect &&
        (!after_end_first->isVisited))
    {
        bool end_forward = decideNextOrderForHole(
            line_vertices[along_end_index], isX,
            isBigger); // !(decideNextOrderForOut(line_vertices[along_end_index], isX, isBigger));
        VNode * hole_next = end_forward ? line_vertices[along_end_index]->next
                                        : line_vertices[along_end_index]->prev;
        if (hole_next->isVisited)
        {
            first_end = along_end_index;
        }
    }

    //deal with the first part
    if (first_start < first_end)
    {
        collectPolygonsInInterval(isX, isBigger, line_vertices, clip_results, first_start,
                                  first_end, first_forward);
    }
    //end of first part
    /**********************************************************************************************/

    /*********************************************************************************************/
    //deal with the second part
    int second_start = along_end_index + 1;
    bool second_forward =
        false; //decideNextOrderForOut(line_vertices[second_start], isX, isBigger);
    VNode * after_end_second =
        forward ? line_vertices[along_end_index]->next : line_vertices[along_end_index]->prev;
    if (line_vertices[along_end_index]->isOnline && after_end_second->isInDirect &&
        (!after_end_second->isVisited))
    {
        bool end_forward = decideNextOrderForHole(
            line_vertices[along_end_index], isX,
            isBigger); //!(decideNextOrderForOut(line_vertices[along_end_index], isX, isBigger));
        VNode * hole_next = end_forward ? line_vertices[along_end_index]->next
                                        : line_vertices[along_end_index]->prev;
        if (!hole_next->isVisited)
        {
            second_start = along_end_index;
            second_forward = forward;
        }
    }

    if (second_start != along_end_index && second_start < end_index)
    {
        second_forward = decideNextOrderForOut(line_vertices[second_start], isX, isBigger);
    }

    if (second_start < end_index)
    {
        collectSecondPart(line_vertices, second_start, end_index, line_vertices[second_start],
                          end_node, second_forward, vertex_count, clip_results, isX, isBigger,
                          secondHalf);
    }
    //end of second part
    /*********************************************************************************************/
}

void RingClipper::construct_clip_rings(vector<VNode *> & firstHalf, vector<VNode *> & secondHalf,
                                       list<RingClipper *> & clip_results)
{
    if (firstHalf.size() + secondHalf.size() < 3)
        return;
    else
    {
        RingClipper * ring_clip = new RingClipper();
        for (auto & vit : firstHalf)
        {
            ring_clip->add(copy(vit));
        }
        for (vector<VNode *>::reverse_iterator rvit = secondHalf.rbegin();
             rvit != secondHalf.rend(); ++rvit)
        {
            ring_clip->add(copy(*rvit));
        }
        ring_clip->rings.push_back(ring_clip->first);

        Ring ring;
        convert_clip_to_ring(ring_clip->first, ring);
        if (fabs(ring.getArea()) < SMALLAREA) //SMALLNUMBER
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
                               bool isNextOrder, list<RingClipper *> & clip_results, bool isX,
                               bool isBigger)
{
    int vsize = line_vertices.size();
    vector<VNode *> firstHalf;
    vector<VNode *> secondHalf;

    //collect first half: out boundary
    VNode * start_node = line_vertices[start_index];
    VNode * end_node = line_vertices[end_index];
    VNode * node = start_node;
    firstHalf.push_back(start_node);
    start_node->isVisited = true;
    do
    {
        if (isNextOrder)
            node = node->next;
        else
            node = node->prev;

        node->isVisited = true;
        firstHalf.push_back(node);
    } while (node != end_node);

    //filter the starting vertex
    int startId = start_index + 1;
    int endId = end_index - 1;
    bool forward = !isNextOrder;
    VNode * start_neighbor =
        !isNextOrder ? line_vertices[start_index]->next : line_vertices[start_index]->prev;
    if (line_vertices[start_index]->isOnline && (!start_neighbor->isVisited))
    {
        assert(start_neighbor->isInDirect);
        forward = !isNextOrder;
        startId = start_index;
    }
    else
    {
        startId = start_index + 1;
        if (startId < endId) forward = decideNextOrderForOut(line_vertices[startId], isX, isBigger);
    }

    //if(line_vertices[end_index]->isOnline

    if (startId < endId)
    {
        collectSecondPart(line_vertices, startId, endId, line_vertices[startId],
                          line_vertices[endId], forward, vsize, clip_results, isX, isBigger,
                          secondHalf);
    }

    //combine two lists to get a new polygon
    construct_clip_rings(firstHalf, secondHalf, clip_results);
}

void RingClipper::collectPolygonsInInterval(bool isX, bool isBigger,
                                            vector<VNode *> & line_vertices,
                                            list<RingClipper *> & clips, int from, int to,
                                            bool start_forward)
{
    int size = line_vertices.size();
    bool forward = start_forward;

    int startId = from;
    while (startId < to)
    {
        VNode * start_node = line_vertices[startId];
        assert(start_node->next->isInDirect || start_node->prev->isInDirect);

        vector<VNode *> along_vertices;
        VNode * end_node = collectAlongWay(along_vertices, start_node, forward);
        int endId = locateVInVector(line_vertices, end_node, startId, size);

        assert(endId != -1);

        collect_ring(line_vertices, startId, endId, forward, clips, isX,
                     isBigger); //note here or not, may have bugs here

        if (endId >= to)
        {
            assert(endId == to);
            break;
        }

        VNode * end_onward = forward ? line_vertices[endId]->next : line_vertices[endId]->prev;
        if (line_vertices[endId]->isOnline && (!end_onward->isVisited))
        {
            startId = endId;
        }
        else
        {
            startId = endId + 1;
            forward = decideNextOrderForOut(line_vertices[startId], isX, isBigger);
        }
    }
}

void RingClipper::collectPolygons(double line_coordinate, bool isX, bool isBigger,
                                  vector<int> & line_vertex_counts, vector<VNode *> & line_vertices,
                                  list<RingClipper *> & clips)
{
    int size = line_vertices.size();
    if (size <= 1)
    {
        if (checkAllIndirect(line_coordinate, isX, isBigger, first))
        {
            clips.push_back(copy_ring_clipper(this, true));
        }
        //set visited
        //for(list<VNode*>::iterator )
        return;
    }

    //sort the intersection vertices
    InsertionSortV(line_vertices, isX);

    bool isNextO = decideNextOrderForOut(line_vertices[0], isX, isBigger); //need modify here

    collectPolygonsInInterval(isX, isBigger, line_vertices, clips, 0, size - 1, isNextO);
}

void RingClipper::clip(double line_coordinate, bool isX, bool isBigger, list<RingClipper *> & clips)
{
    vector<int> line_vertex_counts;
    vector<VNode *> line_vertices;
    //set intersections for holes
    for (auto ring : rings)
    {
        int line_vertex_count =
            setIntersection(ring, line_coordinate, isX, isBigger, line_vertices);
        line_vertex_counts.push_back(line_vertex_count);
    }

    //collect the out boundary
    collectPolygons(line_coordinate, isX, isBigger, line_vertex_counts, line_vertices, clips);

    /**********************************************/
    if (rings.size() == 0 || line_vertex_counts.size() == 0) return;

    list<VNode *>::iterator pit = rings.begin();
    vector<int>::iterator iit = line_vertex_counts.begin();
    //++pit; ++iit; //discard the out boundary
    for (; pit != rings.end(); ++pit, ++iit)
    {

        if (*iit == 0 && ((*pit)->isVisited == false) && (*pit)->isInDirect)
        {
            if ((*pit)->isScanPoint)
            {
                VNode *prev_vertex = nullptr, *next_vertex = nullptr, *prev_adjacent = nullptr,
                      *next_adjacent = nullptr;
                getFeasiblePreNext(*pit, prev_vertex, next_vertex, prev_adjacent, next_adjacent);
                if (!prev_adjacent->isInDirect || !next_adjacent->isInDirect) continue;
            }

            bool haveAdded = false;
            VNode * new_vertex = *pit;
            Ring nrings(Ring::PIN);
            convert_clip_to_ring(new_vertex, nrings);
            const Point2d & center = nrings.getCenter();
            for (auto & clip_part : clips)
            {
                if (clip_ring_encloses(clip_part->first, center))
                {
                    clip_part->rings.push_back(copy_ring(new_vertex));
                    haveAdded = true;
                    break;
                }
            }
            if (!haveAdded)
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
    LinkedVertex * head = ring.getHead();
    VNode * new_head = new VNode(head->getPos()[0], head->getPos()[1]);
    new_head->next = new_head->prev = new_head;
    LinkedVertex * vertex = head->getNext();
    while (vertex != head)
    {
        VNode * node_copy = new VNode(vertex->getPos()[0], vertex->getPos()[1]);
        //leave out a repeated vertex, such as the closing vertex of a shapefile ring
        if (AlmostEqual(node_copy, new_head->prev) || AlmostEqual(node_copy, new_head))
        {
            delete node_copy;
            vertex = vertex->getNext();
            continue;
        }
        new_head->prev->next = node_copy;
        node_copy->prev = new_head->prev;
        node_copy->next = new_head;
        new_head->prev = node_copy;

        vertex = vertex->getNext();
    }

    return new_head;
}

//convert the ring to RingClipper
RingClipper * convert_ring_to_clip(const Ring & ring)
{
    RingClipper * ring_clip = new RingClipper();
    LinkedVertex * head = ring.getHead();
    LinkedVertex * vertex = head;
    do
    {
        VNode * node = new VNode(vertex->getPos()[0], vertex->getPos()[1]);
        ring_clip->add(node);

        vertex = vertex->getNext();
    } while (vertex != head);
    return ring_clip;
}

//convert
void convert_clip_to_ring(VNode * first_vertex, Ring & ring)
{
    ring.beginPoly();
    VNode * vertex = first_vertex;
    do
    {
        ring.addVertex(vertex->x, vertex->y, true);
        vertex = vertex->next;
    } while (vertex != first_vertex);

    ring.endPoly(false);
}

//convert
RingClipper * convertCPolygonToClip(const GIS_polygon & polygon)
{
    RingClipper * ring_clip = new RingClipper();
    for (const auto & part : polygon)
    {
        if (part.getType() == Ring::POUT)
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
void convertClipToPolygon(RingClipper * ring_clip, GIS_polygon * polygon)
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

RingClipper * convertSimplePolygonToClip(const c_polygon & polygon)
{
    RingClipper * ring_clip = new RingClipper();
    for (const auto & part : polygon)
    {
        if (part.getType() == Ring::POUT)
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

VNode * convertBoxToClip(Point2d lower_left, Point2d lower_right, Point2d upper_right,
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

double getOutBdArea(RingClipper * ring_clip)
{
    Ring outer_ring(Ring::POUT);
    convert_clip_to_ring(ring_clip->first, outer_ring);

    double area = outer_ring.getArea();
    outer_ring.destroy();
    return area;
}

//make sure the min and max are in correct order
double intersectTwoBox(Point2d lower_left, Point2d lower_right, Point2d upper_right,
                       Point2d upper_left, double xmin2, double xmax2, double ymin2, double ymax2)
{
    double area(0.0);

    //build a BoxClip
    RingClipper ring_clip;
    VNode * head = convertBoxToClip(lower_left, lower_right, upper_right, upper_left);
    ring_clip.first = head;
    ring_clip.rings.push_back(head);

    list<RingClipper *> clips;
    ring_clip.clip(xmin2, xmax2, ymin2, ymax2, clips);

    for (auto & clip_part : clips)
    {

        area += getOutBdArea(clip_part);

        clip_part->destroy();
        delete clip_part;
    }

    ring_clip.destroy();
    return area;
}

}
