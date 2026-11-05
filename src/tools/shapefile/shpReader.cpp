
#include "shpReader.h"
#include "GridgenError.h"
#include "Warnings.h"
#include "intersection.h"

#include <algorithm>
#include <cassert>
using namespace std;

namespace cusg
{

using uint = unsigned int;

void ShapeReader::destroy()
{
    for (auto & ply : m_ply_list) ply.destroy();
    for (auto & arc : m_arc_list) arc.destroy();
    m_ply_list.clear();
    m_arc_list.clear();
    m_pt_list.clear();

    m_intFlds.clear();
    m_dbFlds.clear();
    m_strFlds.clear();
    m_logFlds.clear();
}

//
//read field attributes from DBF file according to id
bool ShapeReader::readAttributes(DBFHandle & dbfHandle, int iShape)
{
    for (auto & sifld : m_intFlds)
    {
        sifld.vals.insert(make_pair(iShape, DBFReadIntegerAttribute(dbfHandle, iShape, sifld.id)));
    }
    for (auto & sifld : m_dbFlds)
    {
        sifld.vals.insert(make_pair(iShape, DBFReadDoubleAttribute(dbfHandle, iShape, sifld.id)));
    }
    for (auto & sifld : m_strFlds)
    {
        sifld.vals.insert(
            make_pair(iShape, string(DBFReadStringAttribute(dbfHandle, iShape, sifld.id))));
    }
    for (auto & sifld : m_logFlds)
    {
        sifld.vals.insert(
            make_pair(iShape, string(DBFReadLogicalAttribute(dbfHandle, iShape, sifld.id))));
    }
    return true;
}

//
// read from files (both shp and db files)
//
bool ShapeReader::read(const string & filename)
{

    if (already_read) return true;

    //start to read
    SHPHandle hSHP;
    hSHP = SHPOpen(filename.c_str(), "rb");

    //the caller reports a shapefile that cannot be opened
    if (hSHP == nullptr) return false;
    m_filename = filename;

    /************************************************************/
    DBFHandle hDBF = DBFOpen(filename.c_str(), "rb");

    //check
    if (hDBF == nullptr)
    {
        SHPClose(hSHP);
        return false;
    }

    ////check more

    //used for reading from db
    int nWidth, nDecimals;

    //position of these fields in database...
    for (int i = 0; i < DBFGetFieldCount(hDBF); i++)
    {
        char szTitle[12];
        DBFFieldType eType = DBFGetFieldInfo(hDBF, i, szTitle, &nWidth, &nDecimals);
        if (eType != FTInteger && eType != FTDouble && eType != FTString && eType != FTLogical)
            continue; //not good...

        //collect the field name and type info
        if (eType == FTInteger)
        {
            ShapeIntField sintfld;
            sintfld.id = i;
            sintfld.name = szTitle;
            m_intFlds.push_back(sintfld);
        }
        else if (eType == FTDouble)
        {
            ShapeDoubleField sdbfld;
            sdbfld.id = i;
            sdbfld.name = szTitle;
            m_dbFlds.push_back(sdbfld);
        }
        else if (eType == FTString)
        {
            ShapeStringField sstrfld;
            sstrfld.id = i;
            sstrfld.name = szTitle;
            m_strFlds.push_back(sstrfld);
        }
        else if (eType == FTLogical)
        {
            ShapeLogField slogfld;
            slogfld.id = i;
            slogfld.name = szTitle;
            m_logFlds.push_back(slogfld);
        }
    }

    ////check if we can get the field ids for the information we wanted

    /************************************************************/

    //Get info
    int nShapeType, nEntities;
    double adfMinBound[4], adfMaxBound[4];
    SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound);

    //reserve space
    m_ply_list.reserve(nEntities);

    //	Skim over the list of shapes, printing all the vertices
    int nunsupported = 0;
    for (uint i = 0; i < (uint)nEntities; i++)
    {

        //
        SHPObject * psShape = SHPReadObject(hSHP, i);

        if (psShape == nullptr)
        {
            SHPClose(hSHP);
            DBFClose(hDBF);
            throw GridgenError("shape " + to_string(i + 1) + " in shapefile " + filename +
                               " cannot be read");
        }

        //some format check..
        if (psShape->nParts > 0 && psShape->panPartStart[0] != 0)
        {
            cerr << "panPartStart[0] = " << psShape->panPartStart[0] << ", not zero as expected.\n";
        }

        switch (psShape->nSHPType)
        {
            case SHPT_POLYGON:
            case SHPT_POLYGONZ:
            case SHPT_POLYGONM:
            {
                readPly(psShape, i);
                break;
            }
            case SHPT_ARC:
            case SHPT_ARCZ:
            case SHPT_ARCM:
            {
                readArc(psShape, i);
                break;
            }
            case SHPT_POINT:
            case SHPT_POINTM:
            case SHPT_POINTZ:
            case SHPT_MULTIPOINT:
            case SHPT_MULTIPOINTM:
            case SHPT_MULTIPOINTZ:
            {
                readPts(psShape, i);
                break;
            }
            case SHPT_MULTIPATCH:
            {
                readMultiPatch(psShape, i);
                break;
            }
            default:
            {
                nunsupported++;
                break;
            }
        }

        //read the field values
        readAttributes(hDBF, i);

        SHPDestroyObject(psShape);
    }

    SHPClose(hSHP);
    DBFClose(hDBF);

    if (nunsupported > 0)
        warn(to_string(nunsupported) + " of the " + to_string(nEntities) + " shapes in shapefile " +
             filename + " have no geometry or a type that is not supported, and are ignored.");

    //done readh shp file

    //
    //reading from database now.
    //

    already_read = true;
    return true;
}

//whether ring a is inside ring b, which it does not cross; most of the
//vertices of a decide, since a vertex of a may touch b
static bool ringInRing(const vector<Point2d> & a, const vector<Point2d> & b)
{
    int vote = 0;
    for (const Point2d & p : a)
    {
        bool inside = false;
        for (size_t k = 0; k < b.size(); k++)
        {
            if (CrossesRay(b[k].get(), b[(k + 1) % b.size()].get(), p.get())) inside = !inside;
        }
        vote += inside ? 1 : -1;
    }
    return vote > 0;
}

void ShapeReader::readPly(SHPObject * psShape, int iShape)
{
    //the vertices of each part
    vector<vector<Point2d>> rings;
    int iPart = 0;
    for (int j = 0; j < psShape->nVertices; j++)
    {
        if (rings.empty() || (iPart < psShape->nParts && psShape->panPartStart[iPart] == j))
        {
            rings.emplace_back();
            iPart++;
        }

        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        rings.back().emplace_back(x, y);

        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;
    }

    //a ring inside an odd number of the others is a hole of the smallest ring
    //containing it, and any other ring is an outer boundary, whatever the
    //order and orientation of the rings
    size_t nrings = rings.size();
    vector<int> depth(nrings, 0), parent(nrings, -1);
    vector<double> extent(nrings);
    for (size_t i = 0; i < nrings; i++)
    {
        double xmin = rings[i][0][0], xmax = xmin, ymin = rings[i][0][1], ymax = ymin;
        for (const Point2d & p : rings[i])
        {
            xmin = min(xmin, p[0]);
            xmax = max(xmax, p[0]);
            ymin = min(ymin, p[1]);
            ymax = max(ymax, p[1]);
        }
        extent[i] = (xmax - xmin) * (ymax - ymin);
    }
    for (size_t i = 0; i < nrings; i++)
    {
        for (size_t k = 0; k < nrings; k++)
        {
            if (k == i || !ringInRing(rings[i], rings[k])) continue;
            depth[i]++;
            if (parent[i] < 0 || extent[k] < extent[parent[i]]) parent[i] = static_cast<int>(k);
        }
    }

    //the polygon each outer ring starts
    vector<int> polygon_of(nrings, -1);
    for (int hole = 0; hole <= 1; hole++)
    {
        for (size_t i = 0; i < nrings; i++)
        {
            if (depth[i] % 2 != hole) continue;
            c_ply ply(hole ? c_ply::PIN : c_ply::POUT);
            ply.beginPoly();
            for (const Point2d & p : rings[i]) ply.addVertex(p[0], p[1]);
            ply.endPoly();
            if (hole)
                m_ply_list[polygon_of[parent[i]]].push_back(ply);
            else
            {
                polygon_of[i] = static_cast<int>(m_ply_list.size());
                add_ply(ply, iShape);
            }
        }
    }
}
void ShapeReader::readMultiPatch(SHPObject * psShape, int iShape)
{
    //outer and first rings start a polygon, and inner rings are its holes; a
    //plain ring is a hole after another ring
    bool started = false;
    for (int i = 0; i < psShape->nParts; i++)
    {
        int stId = psShape->panPartStart[i];
        int endId = psShape->nVertices - 1;
        if (i < psShape->nParts - 1)
        {
            endId = psShape->panPartStart[i + 1] - 1;
        }
        int type = psShape->panPartType[i];
        if (type != SHPP_OUTERRING && type != SHPP_INNERRING && type != SHPP_FIRSTRING &&
            type != SHPP_RING)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(iShape + 1) +
                 " in shapefile " + m_filename + " is a " + SHPPartTypeName(type) +
                 ", which is not supported, and is skipped.");
            continue;
        }

        if (type == SHPP_INNERRING && !started)
        {
            warn("Part " + to_string(i + 1) + " of multipatch shape " + to_string(iShape + 1) +
                 " in shapefile " + m_filename +
                 " is an inner ring that follows no outer ring, and is skipped.");
            continue;
        }

        bool outer = !started || type == SHPP_OUTERRING || type == SHPP_FIRSTRING;
        c_ply ply(outer ? c_ply::POUT : c_ply::PIN);
        ply.beginPoly();
        for (int j = stId; j <= endId; j++)
        {
            //Points may duplicate....
            double & x = psShape->padfX[j];
            double & y = psShape->padfY[j];
            ply.addVertex(x, y);

            if (x < bbox[0]) bbox[0] = x;
            if (x > bbox[1]) bbox[1] = x;
            if (y < bbox[2]) bbox[2] = y;
            if (y > bbox[3]) bbox[3] = y;
        }
        ply.endPoly();
        add_ply(ply, iShape);
        started = true;
    }
}

void ShapeReader::readArc(SHPObject * psShape, int iShape)
{
    //read each part
    int iPart = 1;

    //
    GIS_plyline arc;
    arc.beginPoly();

    for (int j = 0; j < psShape->nVertices; j++)
    {
        string pszPartType;

        if (j == 0 && psShape->nParts > 0) pszPartType = SHPPartTypeName(psShape->panPartType[0]);

        //this defines a new loop...
        if (iPart < psShape->nParts && psShape->panPartStart[iPart] == j)
        {
            pszPartType = SHPPartTypeName(psShape->panPartType[iPart]);
            iPart++;
            arc.endPoly();
            arc.m_id = iShape;
            arc.computeGeodesicDist();
            m_arc_list.push_back(arc);

            //create a new ply
            arc = GIS_plyline();
            arc.beginPoly();
        }

        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        arc.addVertex(x, y);

        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;
    }

    arc.endPoly();
    arc.m_id = iShape;
    arc.computeGeodesicDist();
    m_arc_list.push_back(arc);
}

void ShapeReader::readPts(SHPObject * psShape, int iShape)
{
    m_pt_list.reserve(psShape->nVertices);
    for (int j = 0; j < psShape->nVertices; j++)
    {
        //Points may duplicate....
        double & x = psShape->padfX[j];
        double & y = psShape->padfY[j];
        if (x < bbox[0]) bbox[0] = x;
        if (x > bbox[1]) bbox[1] = x;
        if (y < bbox[2]) bbox[2] = y;
        if (y > bbox[3]) bbox[3] = y;

        GIS_Point2d pt(x, y);
        pt.m_id = iShape;
        m_pt_list.push_back(pt);
    }
}

void ShapeReader::add_ply(c_ply & ply, int iShape)
{
    if (ply.getType() == c_ply::POUT)
    {
        //create a new polygon
        GIS_polygon polygon;
        polygon.m_id = iShape;
        polygon.push_back(ply);
        m_ply_list.push_back(polygon);
    }
    else
    { //c_ply::PIN
        //add to the last polygon
        m_ply_list.back().push_back(ply);
    }
}

} //end namespace cusg
