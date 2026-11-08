
#include "vert_pass.h"

void determine_area(VPassBox & vertical_pass)
{
    double area = (double)(FLT_MAX);
    double cell_area = fabs(vertical_pass.source->dx * vertical_pass.source->dy);
    if (cell_area < area) area = cell_area;

    cell_area = fabs(vertical_pass.destination->dx * vertical_pass.destination->dy);
    if (cell_area < area) area = cell_area;
    for (auto box : vertical_pass.pass_through)
    {
        cell_area = fabs(box->dx * box->dy);
        if (cell_area < area) area = cell_area;
    }
    vertical_pass.area = area;
}
void gen_one_vp(list<Box *> & column, vector<VPassBox> & passes)
{
    if (column.size() >= 3)
    {
        VPassBox vertical_pass;
        vertical_pass.source = column.front();
        column.pop_front();

        vertical_pass.destination = column.back();
        column.pop_back();

        vertical_pass.pass_through.insert(vertical_pass.pass_through.end(), column.begin(),
                                          column.end());
        determine_area(vertical_pass);

        passes.push_back(vertical_pass);
    }
}
//exploration in two vertical directions
void vertical_explore(Box * box, list<Box *> & column, vector<VPassBox> & passes,
                      bool top_to_down = true)
{
    list<Box *> adjacent_boxes;

    if (top_to_down && box->p_down != nullptr)
        box->p_down->get_leaves(adjacent_boxes);
    else if ((!top_to_down) && (box->p_up != nullptr))
        box->p_up->get_leaves(adjacent_boxes);

    for (auto adjacent : adjacent_boxes)
    {
        assert(adjacent != box);
        if (adjacent == nullptr || !(adjacent->is_leaf)) continue;
        if (adjacent->active)
        {
            if (column.size() >= 2)
            {
                list<Box *> column_with_adjacent(column.begin(), column.end());
                column_with_adjacent.push_back(adjacent);
                gen_one_vp(column_with_adjacent, passes);
            }
        }
        else
        {
            list<Box *> column_with_adjacent(column.begin(), column.end());
            column_with_adjacent.push_back(adjacent);
            vertical_explore(adjacent, column_with_adjacent, passes, top_to_down);
        }
    }
}

void box_vp_connection(Box * b, vector<VPassBox> & passes, bool top_to_down)
{
    if (b == nullptr || (!(b->is_leaf)) || (b->active == false)) return;

    list<Box *> column;
    column.push_back(b);
    vertical_explore(b, column, passes, top_to_down);
}

void get_vertical_passes(Box * b, vector<VPTPair> & pairs, bool top_to_down)
{
    if (b == nullptr || (!(b->is_leaf)) || (b->active == false)) return;

    vector<VPassBox> passes;
    box_vp_connection(b, passes, top_to_down);
    for (vector<VPassBox>::iterator vit = passes.begin(); vit != passes.end(); ++vit)
    {
        Box * source = vit->source;
        Box * destination = vit->destination;

        vector<VPassBox>::iterator jit = vit;
        ++jit;
        for (; jit != passes.end();)
        {
            if (source == jit->source && destination == jit->destination)
            {
                vit->area += jit->area;
                jit = passes.erase(jit);
            }
            else
            {
                ++jit;
            }
        }
    }

    //copy the final results
    for (auto & vertical_pass : passes)
    {
        VPTPair merged(vertical_pass.source, vertical_pass.destination, vertical_pass.area);
        pairs.push_back(merged);
    }

    ////////////////////////////////////////////
}

////this function
