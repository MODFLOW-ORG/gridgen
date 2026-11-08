/*! \file QuadTree3D.cpp

\brief QuadTree 3D grid

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
<a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.

*/

//   $Id: $

#include "QuadTree3D.h"
#include "GridgenError.h"
#include <sstream>
#include "GridIntersection.h"
#include "ascii_grid.h"
#include <algorithm>

#include "vert_pass.h"
#include "RingClipper.h"
using namespace std;

namespace cusg
{

/*!
	Initialize the quadtree grid with a parent grid.  The parent grid should
	be a ModflowGrid object.
	*/

QuadTree3D::QuadTree3D(ModflowGrid * mfgrid)
    : m_box_index_array(nullptr), m_layer_cell_counts(nullptr)
{
    // Store the modflowgrid, too
    m_mfgrid = mfgrid;

    //determine if this is a 2d or 3d grid
    is3d = m_mfgrid->is3d;

    //add the root nodes and their connections.  This creates a node_group
    //object and assigns it to the grid as an attribute.  Nodegroup is
    //an array that contains all of the node objects.
    add_root_nodes();

    //assign a variable that indiactes whether or not box_index_array needs to
    //be rebuilt because nodes have been added or removed.
    m_rebuild_box_index_array = true;
    m_rebuild_cell_array = true;
    m_one_based_node_numbering = true;
    m_box_index_array = get_box_index_array();
    get_cell_array();

    //how z values of top and bot should be determined
    m_bot_layer_z_operator = m_top_layer_z_operator = vector<LayerZValueOperator>(
        m_mfgrid->nlay, Z_VALUE_REPLICATE_OPERATOR); //default is replicate

    vertical_pass_through = false;
}

QuadTree3D::~QuadTree3D(void) = default;

void get_special_point_intersection(QuadTree3D * qtree, const Point2d & pt, int max_level,
                                    int at_layer, list<Box *> & result)
{
    int i, j, x_id, y_id, id;
    ModflowGrid * mf = qtree->get_modflow_grid();
    Index2d id2d = mf->get_index(pt);

    if (id2d[0] == INT_MAX || id2d[1] == INT_MAX || id2d[0] == INT_MIN || id2d[1] == INT_MIN)
        return;

    for (i = -1; i <= 1; i++)
    {
        y_id = id2d[0] + i;
        if (y_id < 0 || y_id >= mf->nrow) continue;

        for (j = -1; j <= 1; j++)
        {
            x_id = id2d[1] + j;
            if (x_id < 0 || x_id >= mf->ncol) continue;

            id = mf->get_box_index(at_layer, y_id, x_id);
            Box * box = qtree->node_group[id];
            if (box->in(pt[0], pt[1]))
            {
                box = box->find(pt[0], pt[1]);
                assert(box);
                result.push_back(box);
            }
        }
    }
}

bool QuadTree3D::refine(Point2d & pt, int max_level, int at_layer)
{
    list<Box *> boxes;
    get_special_point_intersection(this, pt, max_level, at_layer, boxes);
    bool intersects = !boxes.empty();

    while (!boxes.empty())
    {
        Box * b = boxes.front();
        boxes.pop_front();

        if (b->depth >= max_level) continue; //too deep
        bool r = b->split(0);
        assert(r); //make sure that split is successful

        for (auto & c : b->p_children)
        {
            node_group.add(c);
            if (c->in(pt[0], pt[1])) boxes.push_back(c);
        }

        this->m_rebuild_box_index_array = true;
    }

    return intersects;
}

/**************************************************************/

//			//post-check to see if the point is on the corner or on the boundary

//given the  intersection precondition provided by GridIntersection
bool is_really_intersect(Box * b, const LineSeg2d & line)
{
    const double * start = line.pt[0].get();
    const double * end = line.pt[1].get();
    double p[2];
    //create each edge of b
    Point2d corners[4];
    corners[0].set(b->x - b->dx / 2, b->y + b->dy / 2); //upper left corner
    corners[1].set(b->x + b->dx / 2, b->y + b->dy / 2); //upper right corner
    corners[2].set(b->x + b->dx / 2, b->y - b->dy / 2); //lower right corner
    corners[3].set(b->x - b->dx / 2, b->y - b->dy / 2); //lower left corner

    bool b_intersect = false;
    int id1, id2;
    //if collinear
    for (int i = 0; i < 4; i++)
    {
        id1 = i;
        id2 = (i + 1) % 4;
        char code = seg_seg_int(corners[id1].get(), corners[id2].get(), start, end, p);
        if (code == '0')
            continue;
        else if (code == '1')
        {
            b_intersect = true;
            break;
        }
        else if (code == 'v')
        {
            //what vertex is it
            if (between_strict(start, end, corners[id1].get()) ||
                between_strict(start, end, corners[id2].get()))
            {
                b_intersect = true;
                break;
            }
        }
        else if (code == 'e')
        {
            //if any strictly between/inside
            if (between_strict(start, end, corners[id1].get()) ||
                between_strict(start, end, corners[id2].get()) ||
                between_strict(corners[id1].get(), corners[id2].get(), start) ||
                between_strict(corners[id1].get(), corners[id2].get(), end))
            {
                b_intersect = true;
                break;
            }
        }
    }

    if (b_intersect) return true;

    //last chance, both point are in box
    if (b->in(start[0], start[1]) && b->in(end[0], end[1])) return true;

    return false;
}
/***********************************************************************************************************/
//grid and line seg intersection for Quadtree3D refinement
bool special_line_intersect(QuadTree3D * qtree, const LineSeg2d & line, int max_level, int at_layer,
                            vector<Box *> & results)
{
    GridIntersection gi(max_level, at_layer);
    LineSeg2d clipped_line;
    bool clipped = gi.clip_line(qtree, line, clipped_line);
    if (clipped == false)
        return true; //line segment is entirely outside, so there is no intersection
    if (get_sq_length(clipped_line) < 1e-10) return true;

    gi.get_intersections().clear();

    //start from the ends of the clipped segment, since an end outside the grid is in no box
    list<Box *> start_boxes, end_boxes;
    get_special_point_intersection(qtree, clipped_line.pt[0], max_level, at_layer, start_boxes);
    get_special_point_intersection(qtree, clipped_line.pt[1], max_level, at_layer, end_boxes);

    //the boxes at either end, then their neighbors the line crosses
    set<Box *> endpoint_boxes;
    endpoint_boxes.insert(start_boxes.begin(), start_boxes.end());
    endpoint_boxes.insert(end_boxes.begin(), end_boxes.end());

    //a box holding an end point on its boundary may only touch the line
    for (auto b : endpoint_boxes)
    {
        if (is_really_intersect(b, line)) results.push_back(b);
    }

    list<Box *> open(results.begin(), results.end());
    set<Box *> visited(results.begin(), results.end());

    while (open.empty() == false) //not empty
    {
        Box * b = open.front();
        open.pop_front();

        //get the neighbors of the box b
        list<Box *> neighbors[4];

        for (int i = 0; i < 4; i++)
        {
            b->get_neighbors(neighbors[i], i);
        }

        //process neighbors
        for (auto & side : neighbors)
        {
            for (list<Box *>::iterator n = side.begin(); n != side.end(); n++)
            {
                Box * neighbor = *n;

                if (visited.find(neighbor) != visited.end()) continue; //visited

                if (gi.intersect(neighbor, line))
                {
                    if (is_really_intersect(neighbor, line))
                    {
                        results.push_back(neighbor);
                        open.push_back(neighbor);
                    }
                }
                //visited before
                visited.insert(neighbor);
            }
        }

    } //end while

    if (results.empty() == false)
        return true;
    else
        return false;
}

//				if(b->depth>=max_level) continue; //too deep

//	do
//	{

//				if(b->depth>=max_level) continue; //too deep

//				//check if really intersect
//				//todo:

bool QuadTree3D::refine(Polyline & arc, int max_level, int at_layer)
{

    bool refined = true;
    bool intersects = false;

    do
    {

        refined = false;

        vector<Box *> results;

        GridIntersection gi(max_level, at_layer);

        LinkedVertex * ptr = arc.get_head();

        bool line_found = false;
        do
        {
            LinkedVertex * next = ptr->get_next();
            if (next != nullptr)
            {
                LineSeg2d line(ptr->get_pos(), next->get_pos());
                LineSeg2d clipped;
                if (gi.clip_line(this, line, clipped) && get_sq_length(clipped) >= 1e-10)
                    intersects = true;

                if (special_line_intersect(this, line, max_level, at_layer, results))
                    line_found = true;
            }
            ptr = next;

        } while (ptr != nullptr);

        if (line_found)
        {
            for (auto b : results)
            {
                if (b->depth >= max_level) continue; //too deep

                bool was_split = b->split(0);
                if (was_split)
                {
                    for (auto & c : b->p_children)
                    {
                        node_group.add(c);
                    }
                    refined = true;
                    this->m_rebuild_box_index_array = true;
                }
            } //end for i
        }

    } while (refined);

    return intersects;
}

//	do
//	{

bool is_really_intersect(RingClipper * ring_clip, Box * b)
{
    double xmin = min(b->x - b->dx / 2, b->x + b->dx / 2);
    double xmax = max(b->x - b->dx / 2, b->x + b->dx / 2);
    double ymin = min(b->y - b->dy / 2, b->y + b->dy / 2);
    double ymax = max(b->y - b->dy / 2, b->y + b->dy / 2);

    list<RingClipper *> clips;
    ring_clip->clip(xmin, xmax, ymin, ymax, clips);

    double area(0.0);
    for (auto & clip_part : clips)
    {
        GisPolygon clipped_polygon;
        convert_clip_to_polygon(clip_part, &clipped_polygon);
        area += clipped_polygon.get_area();

        clipped_polygon.destroy();
        delete clip_part;
    }

    if (area >= 1e-7)
        return true;
    else
        return false;
}

bool is_really_intersect(Polygon & pgon, Box * b)
{
    RingClipper * ring_clip = convert_simple_polygon_to_clip(pgon);

    bool intersects = is_really_intersect(ring_clip, b);

    ring_clip->destroy();
    delete ring_clip;
    return intersects;
}

//refine cells enclosed by pgon
bool QuadTree3D::refine(Polygon & pgon, int max_level, int at_layer)
{
    bool refined = true;
    bool intersects = false;

    do
    {
        refined = false;
        GridIntersection gi(max_level, at_layer);
        bool r = gi.intersect((QuadTree3D *)this, pgon);
        intersects = intersects || r;

        if (r)
        {
            vector<Box *> & boxes = gi.get_intersections();

            for (auto b : boxes)
            {
                if (b->depth >= max_level) continue; //too deep

                if (!is_really_intersect(pgon, b)) continue; //not really intersect

                if (!b->split(0)) continue; //already split

                for (auto & c : b->p_children) node_group.add(c);

                refined = true;
                this->m_rebuild_box_index_array = true;
            } //end for i (intersected boxes)
        }

    } while (refined);

    return intersects;
}

///return true if this is the final refinement
bool QuadTree3D::smooth_refinement(Box * box, vector<pair<int, Box *>> & box_pq,
                                   int level_tol_horizontal, int level_tol_vertical)
{
    list<pair<int, int>> connection_list = get_node_connections(box);

    bool refine_once = true;

    //check the neighbors
    // go through each connection
    for (auto & connection : connection_list)
    {
        int to_box_index = connection.first;
        int dir = connection.second;
        Box * nei = node_group[to_box_index];

        if (nei == nullptr) continue;
        assert(nei->is_leaf);

        bool refine_nei = false;

        if (abs(dir) == 3) //vertical neighbor
        {
            if (nei->depth < box->depth - level_tol_vertical) refine_nei = true;
            if (nei->depth + 1 < box->depth - level_tol_vertical) refine_once = false;
        }
        else //horizontal neighbor
        {
            if (nei->depth < box->depth - level_tol_horizontal) refine_nei = true;
            if (nei->depth + 1 < box->depth - level_tol_horizontal) refine_once = false;
        }

        //visiting all neighbors
        if (refine_nei)
        {
            bool results = nei->split(0);
            if (results) //the neighbor did split
            {
                for (auto nei_kid : nei->p_children) //enqueue neighbors' kid
                {
                    //add to node group
                    node_group.add(nei_kid);

                    //add to the priority queue
                    box_pq.push_back(make_pair(nei_kid->depth, nei_kid));
                    push_heap(box_pq.begin(), box_pq.end());

                    //flag to rebuild
                    this->m_rebuild_box_index_array = true;
                } //end for k
            }
            else
            {
                throw GridgenError("internal error: a cell cannot be split to smooth the grid");
            }
        }
    } //end for i

    return refine_once;
}

/// Smooth the refinement so no connections exceeds the level tolerance (leveltol).
/// level_tol_horizontal is horizontal level tolerance
/// level_tol_vertical is vertical level tolerance
void QuadTree3D::smooth_refinement(int level_tol_horizontal, int level_tol_vertical)
{
    //get all leaves
    list<Box *> leaves;
    int box_count = node_group.size();
    for (int l = 0; l < box_count; l++)
    {
        Box * box = node_group[l];
        if (box->is_leaf) leaves.push_back(box);
    } //end for l

    //we use a heap here to retrieve the deepest leaves
    vector<pair<int, Box *>> box_pq;
    for (auto box : leaves)
    {
        box_pq.push_back(make_pair(box->depth, box));
        push_heap(box_pq.begin(), box_pq.end());
    }

    //loop until the heap is empty
    while (box_pq.empty() == false)
    {
        Box * box = box_pq.front().second;
        pop_heap(box_pq.begin(), box_pq.end());
        box_pq.pop_back();

        if (box->is_leaf == false) continue; //not a leaf anymore, already split...

        while (smooth_refinement(box, box_pq, level_tol_horizontal, level_tol_vertical) == false);

    } //end while
}

/// Return an array of size (nodes) that contains the area for all cells.
vector<double> QuadTree3D::get_cell_areas()
{
    //log the creation
    //cout<<"Getting cell areas."<<endl;
    int nodesize = nodes();

    vector<double> area(nodesize);

    for (int cell = 0; cell < nodesize; cell++)
    {
        int box_index = m_box_index_array[cell];
        Box * cell_box = node_group[box_index];
        area[cell] = cell_box->dx * cell_box->dy;
    }

    return area;
}

/*!
	\brief Return the four vertices for the specified box_index.
	*/

/*! 
	Return an array of size nodes that contains the box_index for all
	active leaf nodes.  Thus the array contains:

	box_index_array[cell] = box_index

	Note that this array is rebuilt, but only if the
	node_group.rebuild_box_index_array flag is True.  This flag is set to
	true anytime quadtree nodes are added or removed.
	*/
int * QuadTree3D::get_box_index_array()
{
    if (m_box_index_array == nullptr || m_layer_cell_counts == nullptr || m_rebuild_box_index_array)
    {
        node_group.rebuild_box_index_array = true;
        m_box_index_array = node_group.get_box_index_array();
        m_layer_cell_counts = node_group.get_layer_cell_counts();
        m_rebuild_box_index_array = false;
        m_rebuild_cell_array = true;
    }

    return m_box_index_array;
}

///the number of cells in each layer
int * QuadTree3D::get_layer_cell_counts()
{
    if (m_box_index_array == nullptr || m_layer_cell_counts == nullptr || m_rebuild_box_index_array)
    {
        node_group.rebuild_box_index_array = true;
        m_box_index_array = node_group.get_box_index_array();
        m_layer_cell_counts = node_group.get_layer_cell_counts();
        m_rebuild_box_index_array = false;
        m_rebuild_cell_array = true;
    }

    return m_layer_cell_counts;
}

/*!
	\brief Create a cell_array of size node_group with the following:

	a cell is the number of a node in the CSR or COO arrays

	cell_array[box_index] = cell
	*/
int * QuadTree3D::create_cell_array()
{
    int nsize = node_group.size();
    m_cell_array.resize(nsize);
    //fill every byte with that of -INT_MAX, as memset did
    memset(m_cell_array.data(), (int)-INT_MAX, sizeof(int) * nsize);

    int leaf_size = nodes();

    for (int cell = 0; cell < leaf_size; cell++)
    {
        int box_index = m_box_index_array[cell];
        m_cell_array[box_index] = cell;
    }

    return m_cell_array.data();
}

/*! 
	Return an array of len(node_group) that contains the cell for all
	active leaf nodes.  Thus the array contains:

	cell_array[box_index] = cell

	Note that this array is rebuilt, but only if the
	node_group.rebuild_box_index_array flag is True.  This flag is set to
	true anytime quadtree nodes are added or removed.
	*/
int * QuadTree3D::get_cell_array()
{
    if (m_cell_array.empty() || m_rebuild_cell_array)
    {
        create_cell_array();
        m_rebuild_cell_array = false;
    }

    return m_cell_array.data();
}

///get box from id
Box * QuadTree3D::get_cell_box(int cell)
{
    int box_index = m_box_index_array[cell];
    return node_group[box_index];
}

int QuadTree3D::get_cell(int box_index)
{
    int cell = m_cell_array[box_index];
    if (cell == -1)
    {
        throw GridgenError("internal error: cell " + to_string(box_index) + " has no node number");
    }

    return cell;
}

list<pair<int, int>> QuadTree3D::get_node_connections(int cell)
{
    int box_index = m_box_index_array[cell];
    Box * box = node_group[box_index];
    assert(box->is_leaf);
    return get_node_connections(box);
}

list<pair<int, int>> QuadTree3D::get_node_connections(Box * box)
{
    //neighbors from 4 sides and up and down
    list<Box *> nei[6];

    //get each of the 4 neighbors
    for (short i = 0; i < 4; i++) box->get_neighbors(nei[i], i);

    //get all leaves of box->p_up
    if (box->p_up != nullptr)
    {
        box->refine_p_up();
        box->p_up->get_leaves(nei[4]);
    }

    //get all leaves of box->p_down
    if (box->p_down != nullptr)
    {
        box->refine_p_down();
        box->p_down->get_leaves(nei[5]);
    }

    //4 sides and up and down
    list<pair<int, int>> connections;
    for (auto & i : nei[0])
        if (i->active) connections.push_back(make_pair(i->id, 2));
    for (auto & i : nei[3])
        if (i->active) connections.push_back(make_pair(i->id, -1));
    for (auto & i : nei[1])
        if (i->active) connections.push_back(make_pair(i->id, 1));
    for (auto & i : nei[2])
        if (i->active) connections.push_back(make_pair(i->id, -2));
    for (auto & i : nei[4])
        if (i->active) connections.push_back(make_pair(i->id, 3));

    for (auto & i : nei[5])
        if (i->active) connections.push_back(make_pair(i->id, -3));

    return connections;
}

///get box connections from id
///Return a sorted list of connections, where a connection is
///a list [to_box_index, fldir].  fldir is -1, +1, -2, +2 for -x, +x, -y, +y, -z, +z.
list<pair<int, int>> QuadTree3D::get_node_connections_with_pass_through(int cell,
                                                                        vector<double> * area)
{
    int box_index = m_box_index_array[cell];
    Box * box = node_group[box_index];
    assert(box->is_leaf);
    return get_node_connections_with_pass_through(box, area);
}

list<pair<int, int>> QuadTree3D::get_node_connections_with_pass_through(Box * box,
                                                                        vector<double> * area)
{
    //neighbors from 4 sides and up and down
    list<Box *> nei[6];

    //get each of the 4 neighbors
    for (short i = 0; i < 4; i++) box->get_neighbors(nei[i], i);

    //get all leaves of box->p_up
    if (box->p_up != nullptr)
    {
        box->refine_p_up();
        box->p_up->get_leaves(nei[4]);
    }

    //get all leaves of box->p_down
    if (box->p_down != nullptr)
    {
        box->refine_p_down();
        box->p_down->get_leaves(nei[5]);
    }

    //4 sides and up and down
    list<pair<int, int>> connections;
    for (auto & i : nei[0])
        if (i->active) connections.push_back(make_pair(i->id, 2));
    for (auto & i : nei[3])
        if (i->active) connections.push_back(make_pair(i->id, -1));
    for (auto & i : nei[1])
        if (i->active) connections.push_back(make_pair(i->id, 1));
    for (auto & i : nei[2])
        if (i->active) connections.push_back(make_pair(i->id, -2));
    for (auto & i : nei[4])
        if (i->active) connections.push_back(make_pair(i->id, 3));

    //vertical pass here
    if (area != nullptr)
    {
        int csize = connections.size();
        for (int i = 0; i < csize; i++) area->push_back(-1.0f);
    }
    //up direction
    vector<VPTPair> passes_up;
    get_vertical_passes(box, passes_up, false);
    for (vector<VPTPair>::iterator vpit = passes_up.begin(); vpit != passes_up.end(); ++vpit)
    {
        assert(vpit->source == box);
        connections.push_back(make_pair(vpit->destination->id, 3));
        if (area != nullptr) area->push_back(vpit->area);
    }

    for (auto & i : nei[5])
        if (i->active) connections.push_back(make_pair(i->id, -3));

    if (area != nullptr)
    {
        int add_size = connections.size() - area->size();
        for (int i = 0; i < add_size; i++) area->push_back(-1.0f);
    }
    //add the down vertical pass
    vector<VPTPair> passes_down;
    get_vertical_passes(box, passes_down, true);
    for (vector<VPTPair>::iterator vpit = passes_down.begin(); vpit != passes_down.end(); ++vpit)
    {
        assert(vpit->source == box);
        connections.push_back(make_pair(vpit->destination->id, -3));
        if (area != nullptr) area->push_back(vpit->area);
    }

    return connections;
}

/*! \brief Obtain CSRData

	Return a UsgCsrData object, which contains all of the necessary
	information for creating an unstructured grid model.

	If top[nodes] and bot[nodes] are provided they are used in the
	calculation of fahl.

	*/
CSRData QuadTree3D::get_usg_csr_data()
{
    cout << "- Compute the cell connections" << endl;

    if (this->m_rebuild_box_index_array)
    {
        this->get_box_index_array();
    }

    if (this->m_rebuild_cell_array)
    {
        this->get_cell_array();
    }

    //this is the number of leaves of the quadtree
    int cell_count = this->nodes();

    //go through each cell and count diagonal position and connections
    int nja = 0;
    for (int cell = 0; cell < cell_count; cell++)
    {
        nja++;

        vector<double> areas;
        int csize = 0;

        if (vertical_pass_through)
            csize = get_node_connections_with_pass_through(cell, &areas).size();
        else
            csize = get_node_connections(cell).size();
        nja += csize;
    }

    //allocate the CSR arrays.
    //cout<<"Allocating ia of size "<<cell_count+1<<" and ja of size "<<nja<<endl;

    vector<int> ia(cell_count + 1);
    vector<int> ja(nja);
    vector<double> area = get_cell_areas();
    vector<double> cl1(nja);
    vector<double> cl2(nja);
    vector<double> fahl(nja);
    vector<int> fldr(nja);
    vector<int> iac(cell_count); // number of connection +1 for each cell, array of size  nodes

    //Go through each connection and fill the arrays.
    int ja_index = 0;
    int id_offset = (m_one_based_node_numbering) ? 1 : 0;
    for (int cell = 0; cell < cell_count; cell++)
    {
        Box * cell_box = get_cell_box(cell);
        double dxdydzn[3] = {cell_box->dx, cell_box->dy, cell_box->dz};

        //diagonal position
        ia[cell] = ja_index + id_offset;
        ja[ja_index] = cell + id_offset;
        cl1[ja_index] = 0.0;
        cl2[ja_index] = 0.0;
        fahl[ja_index] = 0.0;
        fldr[ja_index] = 0;
        ja_index++;

        vector<double> con_areas;
        list<pair<int, int>> connection_list;
        if (vertical_pass_through)
            connection_list = get_node_connections_with_pass_through(cell, &con_areas);
        else
            connection_list = get_node_connections(cell);

        iac[cell] = connection_list.size() + 1;

        if (sort_connections)
        {
            if (vertical_pass_through)
                sort_areas_by_cell(connection_list, con_areas); //sort the area first
            sort_connections_by_cell(connection_list);
        }

        // go through each connection
        vector<double>::const_iterator fit = con_areas.begin();
        for (auto & connection : connection_list)
        {
            int to_box_index = connection.first;
            int fldir = connection.second;

            Box * to_box = node_group[to_box_index];

            int to_node = get_cell(to_box_index); ///messed up here....

            double to_box_dxdydz[3] = {to_box->dx, to_box->dy, to_box->dz};

            ja[ja_index] = to_node + id_offset;
            fldr[ja_index] = fldir;
            int axis = abs(fldir) - 1;
            cl1[ja_index] = dxdydzn[axis] / 2.0;
            cl2[ja_index] = to_box_dxdydz[axis] / 2.0;
            double dz = 1.0;

            if (is3d)
            {
                dz = 0.5 * (dxdydzn[2] + to_box_dxdydz[2]);
                if (top != nullptr && bot != nullptr)
                {
                    double dz1 = top[cell] - bot[cell];
                    double dz2 = top[to_node] - bot[to_node];
                    dz = 0.5 * (dz1 + dz2);
                }
            }

            if ((!vertical_pass_through) || (*fit) < 0)
            {
                if (axis == 0) //x
                {
                    fahl[ja_index] = min(dxdydzn[1], to_box_dxdydz[1]) * dz;
                }
                else if (axis == 1) //y
                {
                    fahl[ja_index] = min(dxdydzn[0], to_box_dxdydz[0]) * dz;
                }
                else if (axis == 2) //z
                {
                    fahl[ja_index] =
                        min(dxdydzn[0] * dxdydzn[1], to_box_dxdydz[0] * to_box_dxdydz[1]);
                }
            }
            else
            {
                assert(fldir == 3 || fldir == -3);
                fahl[ja_index] = *fit;
            }

            ja_index++;

            if (vertical_pass_through) ++fit;

        } //end for connection

    } //end for cell

    ia[cell_count] = ja_index + id_offset;

    return CSRData(cell_count, nja, std::move(area), std::move(ia), std::move(ja), std::move(cl1),
                   std::move(cl2), std::move(fahl), std::move(fldr), std::move(iac));
}

void QuadTree3D::sort_connections_by_cell(list<pair<int, int>> & connection_list)
{
    //get node numbers
    vector<pair<int, pair<int, int>>> sorted_connection_list;
    for (auto & connection : connection_list)
    {
        int to_node = get_cell(connection.first);
        sorted_connection_list.push_back(make_pair(to_node, connection));
    }
    //sort by node numbers
    sort(sorted_connection_list.begin(), sorted_connection_list.end());
    //put the sorted connection back
    list<pair<int, int>> sorted;
    int size = sorted_connection_list.size();
    for (int i = 0; i < size; i++) sorted.push_back(sorted_connection_list[i].second);
    sorted.swap(connection_list);
}

void QuadTree3D::sort_areas_by_cell(const list<pair<int, int>> & connection_list,
                                    vector<double> & area)
{
    using CIT = list<pair<int, int>>::const_iterator;
    vector<pair<int, double>> sorted_area_list;
    vector<double>::iterator fit = area.begin();
    for (CIT connection = connection_list.begin(); connection != connection_list.end();
         ++connection, ++fit)
    {
        int to_node = get_cell(connection->first);
        sorted_area_list.push_back(make_pair(to_node, *fit));
    }
    //sort by node numbers
    sort(sorted_area_list.begin(), sorted_area_list.end());
    //put the sorted connection back
    vector<double> sorted;
    int size = sorted_area_list.size();
    for (int i = 0; i < size; i++) sorted.push_back(sorted_area_list[i].second);
    sorted.swap(area);
}

void QuadTree3D::add_root_nodes()
{
    // Copy the node_group from the base modflow grid
    node_group.share(m_mfgrid->node_group);
    node_group.nlay = m_mfgrid->nlay;
    node_group.grid = this;

    // Copy the X and Y coordinates for the center of all cells
    // from base modflow grid.
    this->x_centers = clone_array(m_mfgrid->x_centers, m_mfgrid->ncol);
    this->y_centers = clone_array(m_mfgrid->y_centers, m_mfgrid->nrow);

    assert(x_centers && y_centers);
}

// Create the top array from the parent grid
// Do not call this function if you just wanted to get the top array
double * QuadTree3D::get_top()
{
    int nodesize = node_group.node_size();
    int mfgrid_csize = m_mfgrid->nrow * m_mfgrid->ncol;
    //filled in full before it replaces the array, so an error leaves the array unchanged
    unique_ptr<double[]> z(new double[nodesize]);

    for (int cell = 0; cell < nodesize; cell++)
    {
        Box * b = get_cell_box(cell);

        if (m_top_layer_z_operator[b->layer] == Z_VALUE_REPLICATE_OPERATOR)
        {
            Index2d i = m_mfgrid->get_index(Point2d(b->x, b->y));
            z[cell] = m_mfgrid->botm(b->layer, i[0], i[1]);
        }
        else //if(m_top_layer_z_operator[b->layer]==Z_VALUE_LINEAR_INTERPOLATION_OPERATOR)
        {
            if (m_top_layer_z_source.find(b->layer) ==
                m_top_layer_z_source.end()) //no name defined, interpolate from the based grid
                z[cell] = bilinear_interpolation(b, &m_mfgrid->bot[mfgrid_csize * b->layer]);
            else //read from ascii file
            {
                bool area_weighted = true;
                if (m_top_area_weighted_interpolate.find(b->layer) !=
                    m_top_area_weighted_interpolate.end())
                    area_weighted = m_top_area_weighted_interpolate[b->layer];

                z[cell] = bilinear_interpolation(b, m_top_layer_z_source[b->layer], area_weighted);
            }
        }
    }

    delete[] this->top;
    this->top = z.release();
    return this->top;
}

// create the bottom array from the parent grid
// Do not call this function if you just wanted to get the bot array
double * QuadTree3D::get_bot()
{
    int nodesize = nodes();
    int mfgrid_csize = m_mfgrid->nrow * m_mfgrid->ncol;
    //filled in full before it replaces the array, so an error leaves the array unchanged
    unique_ptr<double[]> z(new double[nodesize]);

    for (int cell = 0; cell < nodesize; cell++)
    {
        Box * b = get_cell_box(cell);

        if (m_bot_layer_z_operator[b->layer] == Z_VALUE_REPLICATE_OPERATOR)
        {
            Index2d i = m_mfgrid->get_index(Point2d(b->x, b->y));
            z[cell] = m_mfgrid->botm(b->layer + 1, i[0], i[1]);
        }
        else //if(m_bot_layer_z_operator[b->layer]==Z_VALUE_LINEAR_INTERPOLATION_OPERATOR)
        {
            if (m_bot_layer_z_source.find(b->layer) ==
                m_bot_layer_z_source.end()) //no name defined, interpolate from the based grid
                z[cell] = bilinear_interpolation(b, &m_mfgrid->bot[mfgrid_csize * (b->layer + 1)]);
            else
            {
                bool area_weighted = true;
                if (m_bot_area_weighted_interpolate.find(b->layer) !=
                    m_bot_area_weighted_interpolate.end())
                    area_weighted = m_bot_area_weighted_interpolate[b->layer];

                //read from ascii file
                z[cell] = bilinear_interpolation(b, m_bot_layer_z_source[b->layer], area_weighted);
            }
        }
    }

    delete[] this->bot;
    this->bot = z.release();
    return this->bot;
}

///recompute z value and dz value of each cell based on the top and bottom
bool QuadTree3D::set_cell_elevations()
{
    //check if bot and top are created
    if (bot == nullptr || top == nullptr)
    {
        throw GridgenError("internal error: the cell elevations are used before they are set");
    }

    int nodesize = nodes();

    for (int cell = 0; cell < nodesize; cell++)
    {
        Box * b = get_cell_box(cell);
        b->dz = (top[cell] - bot[cell]);
        b->z = (top[cell] + bot[cell]) / 2;
    }
    return true;
}

///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells)
bool QuadTree3D::align_shared_elevations()
{
    //check if bot and top are created
    if (bot == nullptr || top == nullptr)
    {
        throw GridgenError("internal error: the cell elevations are used before they are set");
    }

    int nodesize = nodes();

    //for each leaf node
    for (int cell = 0; cell < nodesize; cell++)
    {
        Box * box = get_cell_box(cell);
        assert(box->is_leaf); //make sure this node is a leaf node...

        if (box->active == false) continue;   //disabled...
        if (box->p_down == nullptr) continue; //no neighbors at lower level
        align_shared_elevations(box);
    } //end for cell

    return true;
}

///align the elevation (z) between the shared surface (bottom of top cell and top of bottom cells) of the given box
void QuadTree3D::align_shared_elevations(Box * top_box)
{
    //get all leaves of box->p_down
    //TODO: we need to consider VP here...

    double bot_of_top_box = bot[top_box->number];

    list<Box *> lower_boxes;
    top_box->refine_p_down();
    top_box->p_down->get_leaves(lower_boxes);

    //align
    for (auto bot_box : lower_boxes)
    {
        if (bot_box->active == false) continue; //disabled...
        double top_of_bot_box = this->top[bot_box->number];
        if (bot_of_top_box < top_of_bot_box)
        {
            this->top[bot_box->number] = bot_of_top_box;
        }
    } //end for j
}

///number the nodes (i.e. leaves) in the grid
void QuadTree3D::number_nodes()
{
    ModflowGrid * mfgrid = get_modflow_grid();
    unsigned int number = 0;
    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(l, r, c);
                Box * box = node_group[bid];
                number_nodes(box, number);
            } //end for c
        } //end for r
    } //end for l
}

void QuadTree3D::number_nodes(Box * box, unsigned int & number)
{
    if (box->is_leaf)
    {
        if (box->active)
            box->number = number++;
        else
            box->number = -1;
    }
    else
    {
        number_nodes(box->p_children[0], number);
        number_nodes(box->p_children[1], number);
        number_nodes(box->p_children[3], number);
        number_nodes(box->p_children[2], number);
    }
}

//refresh the tree to get the new box index and cell arrays
void QuadTree3D::refresh()
{
    m_rebuild_box_index_array = true;
    m_box_index_array = get_box_index_array();
    get_cell_array();
}

//print the information of the leaves
inline void dfs_print(ostream & out, Box * box, unsigned int & id, const string & extra,
                      string & location, bool one_based_numbering)
{
    char quadrant[8];

    if (box->is_leaf)
    {
        if (box->active)
        {
            if (id != box->number + ((one_based_numbering) ? 1 : 0))
                throw GridgenError("internal error: node " +
                                   to_string(box->number + ((one_based_numbering) ? 1 : 0)) +
                                   " is written as node " + to_string(id));
            out << id++ << ", " << extra << " " << location << "\n";
        }
        else
            out << box->number << ", " << extra << " " << location << "\n";
    }
    else
    {
        for (int x = 0; x < 4; x++)
        {
            int i = 0;
            switch (x)
            {
                case 0: i = 0; break;
                case 1: i = 1; break;
                case 2: i = 3; break;
                case 3: i = 2; break;
            }
            //quadrants are numbered 1 to 4 whatever the node numbering
            sprintf(quadrant, "%d", x + 1);
            string child_location = location + quadrant;
            dfs_print(out, box->p_children[i], id, extra, child_location, one_based_numbering);
        } //end for i
    }
}

ostream & operator<<(ostream & out, QuadTree3D & grid)
{
    ModflowGrid * mfgrid = grid.get_modflow_grid();
    unsigned int id_offset = (grid.m_one_based_node_numbering) ? 1 : 0;
    unsigned int id = id_offset; //starting from 1 if it's one based
    char extra[64];

    out << grid.node_group.get_total_leaf_size() << "\n";

    /* 
		   //
		   // the following code is removed on 3/25/2014. 
		   // this is now output to .nod file
		   // in bool DefinitionExporter::save_node_coordinate(const string& prefix, Grid * grid)
		   // in def_exporter.cpp
		   //

		//output number of nodes and number of connections
		//this is the number of active leaves of the quadtree
		int box_size=grid.nodes();

		//go through each cell and count diagonal position and connections
		int nja = 0;
		for( int nid=0;nid<box_size;nid++)  //nid in range(self.nodes) )
		{
			nja++;
			int csize=grid.get_node_connections(nid).size();
			nja += csize;
		}

		out<<box_size<<"\n";
		out<<nja<<"\n";
		//
		*/

    for (int l = 0; l < mfgrid->nlay; l++)
    {
        for (int r = 0; r < mfgrid->nrow; r++)
        {
            for (int c = 0; c < mfgrid->ncol; c++)
            {
                int bid = mfgrid->get_box_index(l, r, c);
                Box * box = grid.node_group[bid];
                //layer, row, and column are one based whatever the node numbering
                sprintf(extra, "(%d,%d,%d)", l + 1, r + 1, c + 1);
                string extra_str = extra;
                string location;
                dfs_print(out, box, id, extra_str, location, grid.m_one_based_node_numbering);
            } //end for c
        } //end for r
    } //end for l

    return out;
}

void QuadTree3D::save_ghost_node(string file_path)
{
    if (m_gns.size() == 0)
    {
        //computed in full before it is kept, so an error leaves no partial list
        vector<GhostNode> ghost_nodes;
        compute_ghost_nodes(this, ghost_nodes);
        m_gns.swap(ghost_nodes);
    }

    write_ghost_node(file_path, m_gns, ((this->m_one_based_node_numbering) ? 1 : 0));
    cout << "- Save data to file: " << file_path << endl;
}

//linear interpolation the value of the center of the box
double QuadTree3D::bilinear_interpolation(Box * box, double * zgrid)
{
    //create a ascii_grid from grid
    ModflowAsciiGrid ag(this->get_modflow_grid(), zgrid, m_mfgrid->ncol, m_mfgrid->nrow);
    Point2d pt(box->x, box->y);
    return bilinear_interpolation(pt, &ag);
}

//stop with an error for a cell that has no valid value in the ASCII grid
static void nodata_error(const Point2d & center, Box * box, const string & arcinfo_filename)
{
    stringstream msg;
    msg << "ASCII grid " << arcinfo_filename << " has no valid value for the layer "
        << box->layer + 1 << " cell centered at (" << center[0] << ", " << center[1] << ")";
    throw GridgenError(msg.str());
}

//value, or stop with an error if it is the nodata value of ag
static double require_valid(double value, AsciiGrid * ag, const Point2d & center, Box * box,
                            const string & arcinfo_filename)
{
    if (!ag->is_valid(value)) nodata_error(center, box, arcinfo_filename);
    return value;
}

//linear interpolation the value of the center of the box
double QuadTree3D::bilinear_interpolation(Box * box, const string & arcinfo_filename,
                                          bool use_area_weighted)
{
    AsciiGrid * ag = m_str2ag[arcinfo_filename].get();

    if (ag == nullptr) //not found, build from filename
    {
        auto aag = make_unique<ArcinfoAsciiGrid>();
        bool br = aag->build(arcinfo_filename);
        if (br == false)
        {
            throw GridgenError("ASCII grid file " + arcinfo_filename +
                               " does not give a positive NCOLS, NROWS, and CELLSIZE");
        }

        ag = aag.get();
        m_str2ag[arcinfo_filename] = std::move(aag);
    }

    //rotate box center if the grid is rotated...
    Point2d center;
    {
        double cx, cy, rotation;
        get_modflow_grid()->get_rotation(cx, cy, rotation);
        double rx, ry;
        box->rotate_ctr(cx, cy, rotation, rx, ry);
        center.set(rx, ry);
    }

    //#if 0
    if (!use_area_weighted)
    {
        //interpolate here
        return require_valid(bilinear_interpolation(center, ag), ag, center, box, arcinfo_filename);
    }
    //#endif
    //
    //#if 1
    else
    {
        //change to use the area-weighted interpolation
        Point2d corners[4];
        corners[0].set(min(box->x - box->dx / 2, box->x + box->dx / 2),
                       min(box->y - box->dy / 2, box->y + box->dy / 2));
        corners[1].set(max(box->x - box->dx / 2, box->x + box->dx / 2),
                       min(box->y - box->dy / 2, box->y + box->dy / 2));
        corners[2].set(max(box->x - box->dx / 2, box->x + box->dx / 2),
                       max(box->y - box->dy / 2, box->y + box->dy / 2));
        corners[3].set(min(box->x - box->dx / 2, box->x + box->dx / 2),
                       max(box->y - box->dy / 2, box->y + box->dy / 2));

        int dminx = ag->getncols() + 1, dmaxx = -1, dminy = ag->getnrows() + 1, dmaxy = -1;
        double cx, cy, rotation;
        get_modflow_grid()->get_rotation(cx, cy, rotation);
        for (auto & corner : corners)
        {
            double tx, ty;
            tx = cx + (corner[0] - cx) * cos(rotation) - (corner[1] - cy) * sin(rotation);
            ty = cy + (corner[1] - cy) * cos(rotation) + (corner[0] - cx) * sin(rotation);
            corner.set(tx, ty);

            Index2d idx = ag->get_index(corner[0], corner[1]);
            if (idx[1] < dminx) dminx = idx[1];
            if (idx[1] > dmaxx) dmaxx = idx[1];
            if (idx[0] < dminy) dminy = idx[0];
            if (idx[0] > dmaxy) dmaxy = idx[0];
        }

        dminx = max(0, dminx);
        dmaxx = min((int)(ag->getncols() - 1), dmaxx);
        dminy = max(0, dminy);
        dmaxy = min((int)(ag->getnrows() - 1), dmaxy);

        //nodata cells are left out of the average
        double weighted_sum = 0.0, total_area = 0.0, nodata_area = 0.0;
        for (int ix = dminx; ix <= dmaxx; ix++)
        {
            for (int iy = dminy; iy <= dmaxy; iy++)
            {
                AsciiGridCell cell = ag->get_cell(iy, ix);
                double area = intersect_two_box(corners[0], corners[1], corners[2], corners[3],
                                                min(cell.x - cell.dx / 2, cell.x + cell.dx / 2),
                                                max(cell.x - cell.dx / 2, cell.x + cell.dx / 2),
                                                min(cell.y - cell.dy / 2, cell.y + cell.dy / 2),
                                                max(cell.y - cell.dy / 2, cell.y + cell.dy / 2));
                double value = ag->get_value(iy, ix);
                if (!ag->is_valid(value))
                {
                    nodata_area += area;
                    continue;
                }
                weighted_sum += area * value;

                total_area += area;
            }
        }

        //a cell that lies only over nodata cells has no value
        if (total_area < 1e-6 && nodata_area >= 1e-6) nodata_error(center, box, arcinfo_filename);

        //if total intersection area is 0, use the bilinear
        //rotate box center if the grid is rotated...

        if (total_area < 1e-6)
        {
            //interpolate here
            return require_valid(bilinear_interpolation(center, ag), ag, center, box,
                                 arcinfo_filename);
        }
        else
            return weighted_sum / total_area;
    }
    //#endif
}

//// GL: bilinear interpolation

//

//linear interpolation the value of the center of the box
double QuadTree3D::bilinear_interpolation(const Point2d & pt, AsciiGrid * ag)
{
    //determine the i and j to interpolate from...

    //get the ag cell that contains pt
    Index2d id = ag->get_index(pt[0], pt[1]);

    AsciiGridCell cell = ag->get_cell(id[0], id[1]);

    //
    if (pt[0] < cell.x && id[1] > 0) id[1]--;
    if (pt[1] > cell.y && id[0] > 0) id[0]--;

    cell = ag->get_cell(id[0], id[1]);

    //bilinear_interpolation
    double x = (pt[0] - cell.x) / cell.dx;
    double y = (cell.y - pt[1]) / cell.dy;

    // clamp a point outside the cell to its edge
    if (x < 0)
        x = 0;
    else if (x > 1)
        x = 1;
    if (y < 0)
        y = 0;
    else if (y > 1)
        y = 1;

    return cusg::bilinear_interpolation(ag, id[0], id[1], x, y);
}

} //end of namespace
