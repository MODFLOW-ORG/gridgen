
#include "vert_pass.h"

void determineArea(VPassBox & vpb)
{
    double area = (double)(FLT_MAX);
    double cell_area = fabs(vpb.bSrc->dx * vpb.bSrc->dy);
    if (cell_area < area) area = cell_area;

    cell_area = fabs(vpb.bDst->dx * vpb.bDst->dy);
    if (cell_area < area) area = cell_area;
    for (auto tb : vpb.bPssList)
    {
        cell_area = fabs(tb->dx * tb->dy);
        if (cell_area < area) area = cell_area;
    }
    vpb.area = area;
}
void genOneVP(list<Box *> & bl, vector<VPassBox> & vpbl)
{
    if (bl.size() >= 3)
    {
        VPassBox vpb;
        vpb.bSrc = bl.front();
        bl.pop_front();

        vpb.bDst = bl.back();
        bl.pop_back();

        vpb.bPssList.insert(vpb.bPssList.end(), bl.begin(), bl.end());
        determineArea(vpb);

        vpbl.push_back(vpb);
    }
}
//exploration in two vertical directions
void verticalExplore(Box * tb, list<Box *> & bl, vector<VPassBox> & vpbl, bool topToDown = true)
{
    list<Box *> nei;

    if (topToDown && tb->pDown != nullptr)
        tb->pDown->getLeaves(nei);
    else if ((!topToDown) && (tb->pUp != nullptr))
        tb->pUp->getLeaves(nei);

    for (auto pb : nei)
    {
        assert(pb != tb);
        if (pb == nullptr || !(pb->isLeaf)) continue;
        if (pb->active)
        {
            if (bl.size() >= 2)
            {
                list<Box *> tlist(bl.begin(), bl.end());
                tlist.push_back(pb);
                genOneVP(tlist, vpbl);
            }
        }
        else
        {
            list<Box *> tlist(bl.begin(), bl.end());
            tlist.push_back(pb);
            verticalExplore(pb, tlist, vpbl, topToDown);
        }
    }
}

void boxVPConnection(Box * b, vector<VPassBox> & vpbl, bool topToDown)
{
    if (b == nullptr || (!(b->isLeaf)) || (b->active == false)) return;

    list<Box *> bl;
    bl.push_back(b);
    verticalExplore(b, bl, vpbl, topToDown);
}

void getVPT(Box * b, vector<VPTPair> & vpts, bool topToDown)
{
    if (b == nullptr || (!(b->isLeaf)) || (b->active == false)) return;

    vector<VPassBox> vpbs;
    boxVPConnection(b, vpbs, topToDown);
    for (vector<VPassBox>::iterator vit = vpbs.begin(); vit != vpbs.end(); ++vit)
    {
        Box * source = vit->bSrc;
        Box * destination = vit->bDst;

        vector<VPassBox>::iterator jit = vit;
        ++jit;
        for (; jit != vpbs.end();)
        {
            if (source == jit->bSrc && destination == jit->bDst)
            {
                vit->area += jit->area;
                jit = vpbs.erase(jit);
            }
            else
            {
                ++jit;
            }
        }
    }

    //coy the final results
    for (auto & vpb : vpbs)
    {
        VPTPair vpt(vpb.bSrc, vpb.bDst, vpb.area);
        vpts.push_back(vpt);
    }

    ////////////////////////////////////////////
}

////this function
