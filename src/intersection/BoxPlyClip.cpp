#include "BoxPlyClip.h"
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
void InsertionSortV(vector<VNode *> & vs, bool isX)
{
    int vsize = vs.size();
    if (vsize < 2) return;

    for (int i = 1; i < vsize; i++)
    {
        VNode * key = vs[i];
        int j = i - 1;
        while (j >= 0 && sortsAfter(vs[j], key, isX))
        {
            vs[j + 1] = vs[j];
            j--;
        }
        vs[j + 1] = key;
    }
}

int locateVInVector(vector<VNode *> & vs, VNode * vv, int startIdx, int totalSize)
{
    for (int i = startIdx; i < totalSize; i++)
    {
        if (vs[i] == vv) return i;
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

VNode * copyOnePly(VNode * head, bool bSetVisited)
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
BoxPlyClip * copyBoxPlyClip(BoxPlyClip * bpc, bool bSetVisited = false)
{
    BoxPlyClip * clip_copy = new BoxPlyClip();

    VNode * new_head = copyOnePly(bpc->first, bSetVisited);
    clip_copy->first = new_head;
    clip_copy->plys.push_back(new_head);

    if (bpc->plys.size() > 1)
    {
        list<VNode *>::iterator vit = bpc->plys.begin();
        ++vit;
        for (; vit != bpc->plys.end(); ++vit)
        {
            VNode * ply_copy = copyOnePly(*vit, bSetVisited);
            clip_copy->plys.push_back(ply_copy);
        }
    }
    return clip_copy;
}
void BoxPlyClip::resetStatus()
{
    //delete the new added
    for (auto added : newAdds)
    {
        VNode * pre1 = added->prev;
        VNode * next1 = added->next;
        pre1->next = next1;
        next1->prev = pre1;
        added->prev = added->next = nullptr;
        delete added;
    }
    newAdds.clear();

    //loop over the influenced vertices
    for (auto vertex : influVs)
    {
        vertex->resetStatus();
    }
    influVs.clear();
}

double getRotAngle(Vector2d & v1, Vector2d & v2)
{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    //have bug here
    double angle = 0;
    Vector2d dir1 = v1.normalize();
    Vector2d dir2 = v2.normalize();
    double dotpro = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]);
    if (fabs(dotpro - 1.0) <= EPSION)
        return 0;
    else if (fabs(dotpro + 1.0) <= EPSION)
        return PI;
    else
    {
        double cross = 0;
        angle = acos(dotpro);
        cross = (dir1[0] * dir2[1]) - (dir2[0] * dir1[1]);
        if (cross < 0) angle = 2 * PI - angle;

        return angle;
    }
}
//test enclosing
bool ClipPlyEnclose(VNode * head, const Point2d & p2d)
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

void getFeasiblePreNext(VNode * vertex, VNode *& vn1, VNode *& vn2, VNode *& vnprev1,
                        VNode *& vnprev2)
{
    vn1 = vertex->prev;
    vn2 = vertex->next;
    while (AlmostEqual(vn1, vertex))
    {
        vn1 = vn1->prev;
    }
    while (AlmostEqual(vn2, vertex))
    {
        vn2 = vn2->next;
    }
    vnprev1 = vn1;
    vnprev2 = vn2;
    //filtering
    while (vn1->isScanPnt)
    {
        vnprev1 = vn1;
        vn1 = vn1->prev;
    }
    while (vn2->isScanPnt)
    {
        vnprev2 = vn2;
        vn2 = vn2->next;
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
        VNode *vn1 = nullptr, *vn2 = nullptr, *vnpre1 = nullptr, *vnpre2 = nullptr;
        getFeasiblePreNext(vertex, vn1, vn2, vnpre1, vnpre2);
        assert(!(vnpre1->isVisited && vnpre2->isVisited));
        //if one of them is not in indirection
        if (!(vn1->isInDirect && vn2->isInDirect))
        {
            if (vn1->isInDirect && (!vnpre1->isVisited)) //previous
                return false;
            else if (vn2->isInDirect && (!vnpre2->isVisited)) //next
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (vn1->isInDirect && vn2->isInDirect)
        {
            if (vnpre1->isVisited || vn1->isVisited) return true;
            if (vnpre2->isVisited || vn2->isVisited) return false;
        }
        if (!(vnpre1->isScanPnt && vnpre2->isScanPnt))
        {
            Vector2d to_prev(vnpre1->x - vertex->x, vnpre1->y - vertex->y);
            Vector2d to_next(vnpre2->x - vertex->x, vnpre2->y - vertex->y);
            double angle = getRotAngle(to_prev, to_next);
            if ((!isX && isBigger) || (isX && !isBigger))
                return (angle <= PI);
            else
                return (angle >= PI);
        }
        else
        {
            Vector2d to_prev(vnpre1->x - vertex->x, vnpre1->y - vertex->y);
            Vector2d to_next(vn2->x - vertex->x, vn2->y - vertex->y);
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
        VNode *vn1 = nullptr, *vn2 = nullptr, *vnpre1 = nullptr, *vnpre2 = nullptr;
        getFeasiblePreNext(vertex, vn1, vn2, vnpre1, vnpre2);
        assert(!(vnpre1->isVisited && vnpre2->isVisited));
        if (!(vn1->isInDirect && vn2->isInDirect))
        {
            if (vn1->isInDirect && (!vnpre1->isVisited))
                return false;
            else if (vn2->isInDirect && (!vnpre2->isVisited))
                return true;
            else
                assert(false);
        }
        //if both are in direction, but one of them has been visited
        if (vn1->isInDirect && vn2->isInDirect)
        {
            if (vnpre1->isVisited || vn1->isVisited) return true;
            if (vnpre2->isVisited || vn2->isVisited) return false;
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

bool checkAllIndirect(double fval, bool isX, bool isBigger, VNode * head)
{

    Point2d p2d = getCenter(head);
    if (isX)
    {
        return (isBigger == (p2d[0] > fval));
    }
    else
    {
        return (isBigger == (p2d[1] > fval));
    }
}
BoxPlyClip::BoxPlyClip()
{
    first = nullptr;
}
BoxPlyClip::~BoxPlyClip()
{
    destroy();
}
void BoxPlyClip::destroy()
{
    for (auto & ply : plys)
    {
        destroy(ply);
    }
    plys.clear();
    first = nullptr;
}
void BoxPlyClip::destroy(VNode * head)
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

void BoxPlyClip::add(VNode * v)
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
void BoxPlyClip::insert(VNode * before, VNode * after, VNode * v)
{

    before->next = v;
    v->prev = before;
    v->next = after;
    after->prev = v;
}
VNode * BoxPlyClip::collectAlongWay(vector<VNode *> & vs, VNode * start, bool forward)
{
    assert(start->isInDirect && start->isScanPnt);
    VNode * previous = start;
    vs.push_back(previous);
    VNode * node = forward ? start->next : start->prev;

    while (AlmostEqual(previous, node) && node->isInDirect)
    {
        vs.push_back(node);
        previous = node;
        node = forward ? node->next : node->prev;
    }
    ///{
    while (node->isInDirect)
    {
        vs.push_back(node);
        if (node->isScanPnt) break;
        node = forward ? node->next : node->prev;
    }

    return vs.back();
}

//true if v is on the line and both its neighbors are on the kept side, off the line
bool touchesKeptSide(VNode * v)
{
    return v->isOnline && v->next->isInDirect && !v->next->isOnline && v->prev->isInDirect &&
           !v->prev->isOnline;
}

//set intersection
int BoxPlyClip::setIntesection(VNode * startV, double fval, bool isX, bool isBigger,
                               vector<VNode *> & intersectV)
{
    double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
    getXYBound(startV, xmin, xmax, ymin, ymax);
    double xwidth = fabs(xmax - xmin);
    double ywidth = fabs(ymax - ymin);
    if (xwidth < SMALLNUMBER) xwidth = ywidth;
    if (ywidth < SMALLNUMBER) ywidth = xwidth;

    double a[2] = {0}, b[2] = {0};
    if (isX)
    {
        a[0] = fval;
        a[1] = ymin - ywidth / 2;
        b[0] = fval;
        b[1] = ymax + ywidth / 2;
    }
    else
    {
        a[0] = xmin - xwidth / 2;
        a[1] = fval;
        b[0] = xmax + xwidth / 2;
        b[1] = fval;
    }

    //move vertices that are on the line within the tolerance onto it, in place, so that an
    //edge along the line is intersected as one
    VNode * snap = startV;
    do
    {
        if (snap->Online(fval, isX))
        {
            if (isX)
                snap->x = fval;
            else
                snap->y = fval;
        }
        snap = snap->next;
    } while (snap != startV);

    //try to intersect
    map<VNode *, bool> shouldRecord;
    vector<VNode *> interCand;
    VNode * stmp = startV;
    do
    {
        double c[2] = {0}, d[2] = {0}, p[2] = {FLT_MIN, FLT_MIN};
        stmp->setDirection(fval, isX, isBigger);
        //intersect
        c[0] = stmp->x;
        c[1] = stmp->y;
        d[0] = stmp->next->x;
        d[1] = stmp->next->y;
        //an edge with both ends on the line lies along it; SegSegInt can miss this by rounding
        char code = 0;
        if ((isX && c[0] == fval && d[0] == fval) || (!isX && c[1] == fval && d[1] == fval))
            code = 'e';
        else
            code = SegSegInt(a, b, c, d, p);

        //set influenced vertices
        if (stmp->isInDirect) influVs.push_back(stmp);
        if (code == 'e')
        {
            stmp->isOnline = stmp->next->isOnline = true;
            stmp->isScanPnt = stmp->next->isScanPnt = true; //
            influVs.push_back(stmp);
            influVs.push_back(stmp->next);
        }
        else if (code == 'v')
        {
            //if(AlmostEqual(c, p))
            if (AlmostSame(c, p))
            {
                stmp->isOnline = true;
                stmp->isScanPnt = true; //
                influVs.push_back(stmp);
            }
            else
            {
                assert(AlmostSame(d, p));
                stmp->next->isOnline = true;
                stmp->next->isScanPnt = true; //
                influVs.push_back(stmp->next);
            }
        }
        else if (code == '1')
        {
            //if(AlmostEqual(c, p))
            if (AlmostSame(c, p))
            {
                stmp->isOnline = true;
                stmp->isScanPnt = true; //
                influVs.push_back(stmp);
            }
            //else if(AlmostEqual(d, p))
            else if (AlmostSame(d, p))
            {
                stmp->next->isOnline = true;
                stmp->next->isScanPnt = true; //
                influVs.push_back(stmp->next);
            }
            else
            {
                VNode * nv = new VNode(p[0], p[1]);
                this->insert(stmp, stmp->next, nv);
                //add to the newAdds array
                newAdds.push_back(nv);

                nv->isInDirect = true;
                nv->isScanPnt = true; //

                shouldRecord.insert(make_pair(nv, true));
                stmp = stmp->next;
            }
        }

        stmp = stmp->next;
    } while (stmp != startV);

    //loop to decide online v
    stmp = startV;
    do
    {
        if (stmp->isOnline)
        {
            stmp->isInDirect = true;
            interCand.push_back(stmp);
            shouldRecord.insert(make_pair(stmp, true));
        }
        stmp = stmp->next;
    } while (stmp != startV);

    vector<VNode *> uselessNodes;
    vector<VNode *> notOnline;
    vector<VNode *> touching;
    for (auto vertex : interCand)
    {
        bool pureInDir1 = (vertex->next->isInDirect && (!vertex->next->isOnline));
        bool pureInDir2 = (vertex->prev->isInDirect && (!vertex->prev->isOnline));
        if (vertex->isOnline && (!(pureInDir1 && pureInDir2)))
        {
            notOnline.push_back(vertex);
        }

        //a vertex that touches the line from the kept side does not cross it, unless the
        //ring also crosses the line at the same place
        if (touchesKeptSide(vertex))
        {
            bool alone = true;
            for (auto & mit : shouldRecord)
            {
                if (mit.first != vertex && AlmostEqual(mit.first, vertex) &&
                    !touchesKeptSide(mit.first))
                    alone = false;
            }
            if (alone)
            {
                touching.push_back(vertex);
                shouldRecord[vertex] = false;
                continue;
            }
        }

        if (pureInDir1 || pureInDir2)
            shouldRecord[vertex] = true;
        else
            shouldRecord[vertex] = false;

        if (!(pureInDir1 || pureInDir2))
        {
            uselessNodes.push_back(vertex);
        }
    }
    //not online
    for (auto vertex : notOnline)
    {
        vertex->isOnline = false;
    }
    //useless nodes
    for (auto vertex : uselessNodes)
    {
        vertex->isOnline = false;
        vertex->isInDirect = false; // vertex->isScanPnt = false;
    }
    //touching nodes are kept as ordinary vertices
    for (auto vertex : touching)
    {
        vertex->isOnline = false;
        vertex->isScanPnt = false;
    }

    //decide the final intersection vertices, in ring order rather than by address so that
    //vertices at the same place on the scanning line are always sorted the same way
    int onNum = 0;
    stmp = startV;
    do
    {
        map<VNode *, bool>::iterator mit = shouldRecord.find(stmp);
        if (mit != shouldRecord.end() && mit->second)
        {
            intersectV.push_back(stmp);
            onNum++;
        }
        stmp = stmp->next;
    } while (stmp != startV);

    return onNum;
}

void BoxPlyClip::collectSecondPart(vector<VNode *> & interVs, int sIdx, int eIdx,
                                   VNode * start_node, VNode * end_node, bool isNextOrder,
                                   int vsize, list<BoxPlyClip *> & bpcResults, bool isX,
                                   bool isBigger, vector<VNode *> & secondHalf)
{
    int startId = sIdx;
    bool forward = isNextOrder;

    //collect these vertices
    vector<VNode *> tmpVec;
    VNode * tsvTM = collectAlongWay(tmpVec, interVs[startId], forward);

    for (auto vertex : tmpVec)
    {
        if (!vertex->isVisited)
        {
            secondHalf.push_back(vertex);
            vertex->isVisited = true;
        }
    }

    int tsvEndIdx = locateVInVector(interVs, tsvTM, startId, vsize);
    assert(tsvEndIdx != -1);

    /**************************************************************/
    //deal with the first part
    int fsIdx = startId + 1;
    bool fnextO = false; //decideNextOrderForOut(interVs[fsIdx], isX, isBigger);
    VNode * start_neighbor = forward ? interVs[startId]->prev : interVs[startId]->next;
    if (interVs[startId]->isOnline && start_neighbor->isInDirect && (!start_neighbor->isVisited))
    {
        fsIdx = startId;
        fnextO = !forward;
    }
    if (fsIdx != startId && fsIdx < tsvEndIdx)
        fnextO = decideNextOrderForOut(interVs[fsIdx], isX, isBigger);

    int feIdx = tsvEndIdx - 1;
    VNode * end_neighbor = forward ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
    if (interVs[tsvEndIdx]->isOnline && end_neighbor->isInDirect && (!end_neighbor->isVisited))
    {
        bool end_forward = decideNextOrderForHole(
            interVs[tsvEndIdx], isX,
            isBigger); // !(decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger));
        VNode * enext = end_forward ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
        if (enext->isVisited)
        {
            feIdx = tsvEndIdx;
        }
    }

    //deal with the first part
    if (fsIdx < feIdx)
    {
        collectPolygonsInInterval(isX, isBigger, interVs, bpcResults, fsIdx, feIdx, fnextO);
    }
    //end of first part
    /**********************************************************************************************/

    /*********************************************************************************************/
    //deal with the second part
    int ssIdx = tsvEndIdx + 1;
    bool sorder = false; //decideNextOrderForOut(interVs[ssIdx], isX, isBigger);
    VNode * vn1 = forward ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
    if (interVs[tsvEndIdx]->isOnline && vn1->isInDirect && (!vn1->isVisited))
    {
        bool end_forward = decideNextOrderForHole(
            interVs[tsvEndIdx], isX,
            isBigger); //!(decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger));
        VNode * snext = end_forward ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
        if (!snext->isVisited)
        {
            ssIdx = tsvEndIdx;
            sorder = forward;
        }
    }

    if (ssIdx != tsvEndIdx && ssIdx < eIdx)
    {
        sorder = decideNextOrderForOut(interVs[ssIdx], isX, isBigger);
    }

    if (ssIdx < eIdx)
    {
        collectSecondPart(interVs, ssIdx, eIdx, interVs[ssIdx], end_node, sorder, vsize, bpcResults,
                          isX, isBigger, secondHalf);
    }
    //end of second part
    /*********************************************************************************************/
}

void BoxPlyClip::constructClipPly(vector<VNode *> & firstHalf, vector<VNode *> & secondHalf,
                                  list<BoxPlyClip *> & bpcResults)
{
    if (firstHalf.size() + secondHalf.size() < 3)
        return;
    else
    {
        BoxPlyClip * bpc = new BoxPlyClip();
        for (auto & vit : firstHalf)
        {
            bpc->add(copy(vit));
        }
        for (vector<VNode *>::reverse_iterator rvit = secondHalf.rbegin();
             rvit != secondHalf.rend(); ++rvit)
        {
            bpc->add(copy(*rvit));
        }
        bpc->plys.push_back(bpc->first);

        c_ply cp;
        convertClipToPly(bpc->first, cp);
        if (fabs(cp.getArea()) < SMALLAREA) //SMALLNUMBER
        {
            cp.destroy();
            delete bpc;
            return;
        }
        else
        {
            cp.destroy();
            bpcResults.push_back(bpc);
        }
    }
}
void BoxPlyClip::collectOnePly(vector<VNode *> & interVs, int sIdx, int eIdx, bool isNextOrder,
                               list<BoxPlyClip *> & bpcResults, bool isX, bool isBigger)
{
    int vsize = interVs.size();
    vector<VNode *> firstHalf;
    vector<VNode *> secondHalf;

    //collect first half: out boundary
    VNode * start_node = interVs[sIdx];
    VNode * end_node = interVs[eIdx];
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
    int startId = sIdx + 1;
    int endId = eIdx - 1;
    bool forward = !isNextOrder;
    VNode * vno = !isNextOrder ? interVs[sIdx]->next : interVs[sIdx]->prev;
    if (interVs[sIdx]->isOnline && (!vno->isVisited))
    {
        assert(vno->isInDirect);
        forward = !isNextOrder;
        startId = sIdx;
    }
    else
    {
        startId = sIdx + 1;
        if (startId < endId) forward = decideNextOrderForOut(interVs[startId], isX, isBigger);
    }

    //if(interVs[eIdx]->isOnline

    if (startId < endId)
    {
        collectSecondPart(interVs, startId, endId, interVs[startId], interVs[endId], forward, vsize,
                          bpcResults, isX, isBigger, secondHalf);
    }

    //combine two lists to get a new polygon
    constructClipPly(firstHalf, secondHalf, bpcResults);
}

void BoxPlyClip::collectPolygonsInInterval(bool isX, bool isBigger, vector<VNode *> & interVs,
                                           list<BoxPlyClip *> & cres, int from, int to, bool ino)
{
    int size = interVs.size();
    bool isNextO = ino;

    int startId = from;
    while (startId < to)
    {
        VNode * start_node = interVs[startId];
        assert(start_node->next->isInDirect || start_node->prev->isInDirect);

        vector<VNode *> tvec;
        VNode * end_node = collectAlongWay(tvec, start_node, isNextO);
        int endId = locateVInVector(interVs, end_node, startId, size);

        assert(endId != -1);

        collectOnePly(interVs, startId, endId, isNextO, cres, isX,
                      isBigger); //note here or not, may have bugs here

        if (endId >= to)
        {
            assert(endId == to);
            break;
        }

        VNode * vn1 = isNextO ? interVs[endId]->next : interVs[endId]->prev;
        if (interVs[endId]->isOnline && (!vn1->isVisited))
        {
            startId = endId;
        }
        else
        {
            startId = endId + 1;
            isNextO = decideNextOrderForOut(interVs[startId], isX, isBigger);
        }
    }
}

void BoxPlyClip::collectPolygons(double fval, bool isX, bool isBigger, vector<int> & voNums,
                                 vector<VNode *> & interVs, list<BoxPlyClip *> & cres)
{
    int size = interVs.size();
    if (size <= 1)
    {
        if (checkAllIndirect(fval, isX, isBigger, first))
        {
            cres.push_back(copyBoxPlyClip(this, true));
        }
        //set visited
        //for(list<VNode*>::iterator )
        return;
    }

    //sort the intersection vertices
    InsertionSortV(interVs, isX);

    bool isNextO = decideNextOrderForOut(interVs[0], isX, isBigger); //need modify here

    collectPolygonsInInterval(isX, isBigger, interVs, cres, 0, size - 1, isNextO);
}

void BoxPlyClip::clip(double fval, bool isX, bool isBigger, list<BoxPlyClip *> & cres)
{
    vector<int> onlineNums;
    vector<VNode *> interVs;
    //set intersections for holes
    for (auto ply : plys)
    {
        int vo = setIntesection(ply, fval, isX, isBigger, interVs);
        onlineNums.push_back(vo);
    }

    //collect the out boundary
    collectPolygons(fval, isX, isBigger, onlineNums, interVs, cres);

    /**********************************************/
    if (plys.size() == 0 || onlineNums.size() == 0) return;

    list<VNode *>::iterator pit = plys.begin();
    vector<int>::iterator iit = onlineNums.begin();
    //++pit; ++iit; //discard the out boundary
    for (; pit != plys.end(); ++pit, ++iit)
    {

        if (*iit == 0 && ((*pit)->isVisited == false) && (*pit)->isInDirect)
        {
            if ((*pit)->isScanPnt)
            {
                VNode *vn1 = nullptr, *vn2 = nullptr, *vnprev1 = nullptr, *vnprev2 = nullptr;
                getFeasiblePreNext(*pit, vn1, vn2, vnprev1, vnprev2);
                if (!vnprev1->isInDirect || !vnprev2->isInDirect) continue;
            }

            bool haveAdded = false;
            VNode * nv = *pit;
            c_ply nply(c_ply::PIN);
            convertClipToPly(nv, nply);
            const Point2d & ctr = nply.getCenter();
            for (auto & cre : cres)
            {
                if (ClipPlyEnclose(cre->first, ctr))
                {
                    cre->plys.push_back(copyOnePly(nv));
                    haveAdded = true;
                    break;
                }
            }
            if (!haveAdded)
            {
                BoxPlyClip * bpc = convertPlyToClip(nply);
                bpc->plys.push_back(bpc->first);
                cres.push_back(bpc);
            }
            nply.destroy(); //clear memory
        }
    }
}
void BoxPlyClip::clip(double xmin, double xmax, double ymin, double ymax, list<BoxPlyClip *> & cres)
{
    //pre-processing
    list<BoxPlyClip *> res1, res2, res3;
    this->clip(xmin, true, true, res1);

    for (auto bpc : res1)
    {
        bpc->clip(xmax, true, false, res2);

        bpc->destroy();
        delete bpc;
    }

    for (auto bpc : res2)
    {
        bpc->clip(ymin, false, true, res3);

        bpc->destroy();
        delete bpc;
    }

    for (auto bpc : res3)
    {
        bpc->clip(ymax, false, false, cres);

        bpc->destroy();
        delete bpc;
    }
}

VNode * BoxPlyClip::copy(const VNode * v)
{
    VNode * nv = new VNode(v->x, v->y);
    return nv;
}

VNode * convertPlyToVN(const c_ply & tply)
{
    ply_vertex * head = tply.getHead();
    VNode * new_head = new VNode(head->getPos()[0], head->getPos()[1]);
    new_head->next = new_head->prev = new_head;
    ply_vertex * vertex = head->getNext();
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

//convert the ply to BoxPlyClip
BoxPlyClip * convertPlyToClip(const c_ply & tply)
{
    BoxPlyClip * bpc = new BoxPlyClip();
    ply_vertex * thead = tply.getHead();
    ply_vertex * vertex = thead;
    do
    {
        VNode * tnode = new VNode(vertex->getPos()[0], vertex->getPos()[1]);
        bpc->add(tnode);

        vertex = vertex->getNext();
    } while (vertex != thead);
    return bpc;
}

//convert
void convertClipToPly(VNode * fv, c_ply & resPly)
{
    resPly.beginPoly();
    VNode * nf = fv;
    do
    {
        resPly.addVertex(nf->x, nf->y, true);
        nf = nf->next;
    } while (nf != fv);

    resPly.endPoly(false);
}

//convert
BoxPlyClip * convertCPolygonToClip(const GIS_polygon & polygon)
{
    BoxPlyClip * bpc = new BoxPlyClip();
    for (const auto & ccp : polygon)
    {
        if (ccp.getType() == c_ply::POUT)
        {
            VNode * nvn = convertPlyToVN(ccp);
            bpc->first = nvn;
            bpc->plys.push_front(nvn);
        }
        else
        {
            VNode * nvn = convertPlyToVN(ccp);
            //a hole with fewer than three distinct vertices has no area
            if (nvn->next == nvn || nvn->next->next == nvn)
                bpc->destroy(nvn);
            else
                bpc->plys.push_back(nvn);
        }
    }
    return bpc;
}
//convert
void convertClipToPolygon(BoxPlyClip * bpc, GIS_polygon * cpolygon)
{
    if (bpc == nullptr || bpc->first == nullptr) return;
    //loop over
    c_ply outply(c_ply::POUT);
    convertClipToPly(bpc->first, outply);
    cpolygon->push_back(outply);

    if (bpc->plys.size() == 1) return;
    list<VNode *>::iterator vit = bpc->plys.begin();
    ++vit;
    for (; vit != bpc->plys.end(); ++vit)
    {
        c_ply tp(c_ply::PIN);
        convertClipToPly(*vit, tp);
        cpolygon->push_back(tp);
    }
}

BoxPlyClip * convertSimplePolygonToClip(const c_polygon & polygon)
{
    BoxPlyClip * bpc = new BoxPlyClip();
    for (const auto & ccp : polygon)
    {
        if (ccp.getType() == c_ply::POUT)
        {
            VNode * nvn = convertPlyToVN(ccp);
            bpc->first = nvn;
            bpc->plys.push_front(nvn);
        }
        else
        {
            VNode * nvn = convertPlyToVN(ccp);
            //a hole with fewer than three distinct vertices has no area
            if (nvn->next == nvn || nvn->next->next == nvn)
                bpc->destroy(nvn);
            else
                bpc->plys.push_back(nvn);
        }
    }
    return bpc;
}

VNode * convertBoxToClip(Point2d ld, Point2d rd, Point2d ru, Point2d lu)
{
    VNode * tn[4];
    tn[0] = new VNode(ld[0], ld[1]);
    tn[1] = new VNode(rd[0], rd[1]);
    tn[2] = new VNode(ru[0], ru[1]);
    tn[3] = new VNode(lu[0], lu[1]);

    VNode * head = tn[0];
    head->next = head->prev = tn[0];

    for (int i = 1; i <= 3; i++)
    {
        head->prev->next = tn[i];
        tn[i]->prev = head->prev;
        tn[i]->next = head;
        head->prev = tn[i];
    }

    return head;
}

double getOutBdArea(BoxPlyClip * bpc)
{
    c_ply outply(c_ply::POUT);
    convertClipToPly(bpc->first, outply);

    double area = outply.getArea();
    outply.destroy();
    return area;
}

//make sure the min and max are in correct order
double intersectTwoBox(Point2d ld, Point2d rd, Point2d ru, Point2d lu, double xmin2, double xmax2,
                       double ymin2, double ymax2)
{
    double area(0.0);

    //build a BoxClip
    BoxPlyClip bpc;
    VNode * head = convertBoxToClip(ld, rd, ru, lu);
    bpc.first = head;
    bpc.plys.push_back(head);

    list<BoxPlyClip *> cres;
    bpc.clip(xmin2, xmax2, ymin2, ymax2, cres);

    for (auto & cre : cres)
    {

        area += getOutBdArea(cre);

        cre->destroy();
        delete cre;
    }

    bpc.destroy();
    return area;
}

}
