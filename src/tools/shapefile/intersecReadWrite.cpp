#include "intersecReadWrite.h"
#include "RoundTrip.h"
#include <sstream>
using namespace std;
using namespace cusg;

const string _IFO_DEL = ",";

// signed area of a ring, positive when it runs counterclockwise
static double ringArea(ply_vertex * head)
{
    double area = 0.0;
    ply_vertex * pv = head;
    do
    {
        const Point2d & a = pv->getPos();
        const Point2d & b = pv->getNext()->getPos();
        area += a[0] * b[1] - b[0] * a[1];
        pv = pv->getNext();
    } while (pv != head);
    return 0.5 * area;
}

bool checkAttributesContain(const char * attribute_name, vector<string> & attributes)
{
    for (auto & str : attributes)
    {
        if (str == attribute_name || str == "all")
            return true;
        else if (str == "none")
            return false;
    }
    return false;
}
void checkShpAttrContain(GridIntersection * grid_intersection, vector<bool> & bIntFld,
                         vector<bool> & bDBFld, vector<bool> & bStrFld, vector<bool> & bLogFld,
                         vector<string> & iName, vector<string> & dName, vector<string> & sName,
                         vector<string> & lName, vector<string> & attributes)
{
    vector<ShapeIntField> & m_intFlds = grid_intersection->m_shp_reader.getIntFields();
    vector<ShapeDoubleField> & m_dbFlds = grid_intersection->m_shp_reader.getDBFields();
    vector<ShapeStringField> & m_strFlds = grid_intersection->m_shp_reader.getStrFields();
    vector<ShapeLogField> & m_logFlds = grid_intersection->m_shp_reader.getLogicalFields();
    for (auto & fld : m_intFlds)
    {
        if (checkAttributesContain(fld.name.c_str(), attributes))
        {
            bIntFld.push_back(true);
            iName.push_back(fld.name);
        }
        else
            bIntFld.push_back(false);
    }
    for (auto & fld : m_dbFlds)
    {
        if (checkAttributesContain(fld.name.c_str(), attributes))
        {
            bDBFld.push_back(true);
            dName.push_back(fld.name);
        }
        else
            bDBFld.push_back(false);
    }
    for (auto & fld : m_strFlds)
    {
        if (checkAttributesContain(fld.name.c_str(), attributes))
        {
            bStrFld.push_back(true);
            sName.push_back(fld.name);
        }
        else
            bStrFld.push_back(false);
    }
    for (auto & fld : m_logFlds)
    {
        if (checkAttributesContain(fld.name.c_str(), attributes))
        {
            bLogFld.push_back(true);
            lName.push_back(fld.name);
        }
        else
            bLogFld.push_back(false);
    }
}
void addShpAttrVal(int iShape, GridIntersection * grid_intersection, vector<bool> & bIntFld,
                   vector<bool> & bDBFld, vector<bool> & bStrFld, vector<bool> & bLogFld,
                   vector<int> & int_values, vector<double> & double_values,
                   vector<string> & string_values, vector<string> & logical_values)
{
    vector<ShapeIntField> & m_intFlds = grid_intersection->m_shp_reader.getIntFields();
    vector<ShapeDoubleField> & m_dbFlds = grid_intersection->m_shp_reader.getDBFields();
    vector<ShapeStringField> & m_strFlds = grid_intersection->m_shp_reader.getStrFields();
    vector<ShapeLogField> & m_logFlds = grid_intersection->m_shp_reader.getLogicalFields();
    vector<bool>::iterator ibit = bIntFld.begin();
    for (vector<ShapeIntField>::iterator sit = m_intFlds.begin(); sit != m_intFlds.end();
         ++sit, ++ibit)
    {
        if (*ibit) int_values.push_back(sit->values[iShape]);
    }
    vector<bool>::iterator idit = bDBFld.begin();
    for (vector<ShapeDoubleField>::iterator sit = m_dbFlds.begin(); sit != m_dbFlds.end();
         ++sit, ++idit)
    {
        if (*idit) double_values.push_back(sit->values[iShape]);
    }
    vector<bool>::iterator sbit = bStrFld.begin();
    for (vector<ShapeStringField>::iterator sit = m_strFlds.begin(); sit != m_strFlds.end();
         ++sit, ++sbit)
    {
        if (*sbit) string_values.push_back(sit->values[iShape]);
    }
    vector<bool>::iterator lbit = bLogFld.begin();
    for (vector<ShapeLogField>::iterator sit = m_logFlds.begin(); sit != m_logFlds.end();
         ++sit, ++lbit)
    {
        if (*lbit) logical_values.push_back(sit->values[iShape]);
    }
}

void addDBHandleFields(DBFHandle & dbfHandle, int intFldNum, vector<string> & intFldNames,
                       int doubleFldNum, vector<string> & dbFldNames, vector<int> & intFldIDs,
                       vector<int> & dbFldIDs)
{
    for (int i = 0; i < intFldNum; i++)
    {
        int tid = DBFAddField(dbfHandle, intFldNames[i].c_str(), FTInteger, 20, 0);
        intFldIDs.push_back(tid);
    }
    for (int i = 0; i < doubleFldNum; i++)
    {
        int tid = DBFAddField(dbfHandle, dbFldNames[i].c_str(), FTDouble, 20, 8);
        dbFldIDs.push_back(tid);
    }
}
void addDBHandleStrFields(DBFHandle & dbfHandle, int strFldNum, vector<string> & strFldNames,
                          vector<int> & strFldIDs)
{
    for (int i = 0; i < strFldNum; i++)
    {
        int tid = DBFAddField(dbfHandle, strFldNames[i].c_str(), FTString, 50, 0);
        strFldIDs.push_back(tid);
    }
}
//write int and double attributes
void writeIntDBVal(DBFHandle & dbfHandle, int entity, vector<int> & intFldIDs,
                   vector<int> & intVals, vector<int> & dbFldIDs, vector<double> & dbVals,
                   vector<int> & strFldIDs, vector<string> & strVals, vector<int> & logFldIDs,
                   vector<string> & logVals, ofstream & ofile, bool isTxt)
{
    int intFldNum = intVals.size();
    int dbFldNum = dbVals.size();
    int strFldNum = strVals.size();
    int logFldNum = logVals.size();
    //write integer fields
    for (int i = 0; i < intFldNum; ++i)
    {
        if (!isTxt)
            DBFWriteIntegerAttribute(dbfHandle, entity, intFldIDs[i], intVals[i]);
        else
            ofile << intVals[i] << _IFO_DEL;
    }
    for (int i = 0; i < dbFldNum; ++i)
    {
        if (!isTxt)
            DBFWriteDoubleAttribute(dbfHandle, entity, dbFldIDs[i], dbVals[i]);
        else
            ofile << round_trip(dbVals[i]) << _IFO_DEL;
    }
    for (int i = 0; i < strFldNum; i++)
    {
        if (!isTxt)
            DBFWriteStringAttribute(dbfHandle, entity, strFldIDs[i], strVals[i].c_str());
        else
            ofile << strVals[i].c_str() << _IFO_DEL;
    }
    for (int i = 0; i < logFldNum; i++)
    {
        if (!isTxt)
            DBFWriteStringAttribute(dbfHandle, entity, logFldIDs[i], logVals[i].c_str());
        else
            ofile << logVals[i].c_str() << _IFO_DEL;
    }
    ofile << endl;
}
void writeIntersection(string file, int intFldNum, vector<vector<int>> & intfldVals,
                       vector<string> & intFldNames, int doubleFldNum,
                       vector<vector<double>> & dbfldVals, vector<string> & dbFldNames,
                       int strFldNum, vector<vector<string>> & strfldVals,
                       vector<string> & strFldNames, int logFldNum,
                       vector<vector<string>> & logfldVals, vector<string> & logFldNames,
                       GridIntersection * grid_intersection, WRITETYPE wtype, bool isTxt)
{
    //create a shape file
    ofstream ofile;
    SHPHandle handler = nullptr;
    if (wtype == WRITE_POINT)
    {
        if (!isTxt)
            handler = SHPCreate(file.c_str(), SHPT_POINT);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }
    else if (wtype == WRITE_POLYLINE)
    {
        if (!isTxt)
            handler = SHPCreate(file.c_str(), SHPT_ARC);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }
    else
    {
        if (!isTxt)
            handler = SHPCreate(file.c_str(), SHPT_POLYGON);
        else
        {
            ofile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }

    //create a .dbf file
    DBFHandle dbfHandle = nullptr;
    if (!isTxt)
    {
        string dbf = file;
        if (file.find('.') != string::npos)
        {
            int dot_index = file.find_last_of('.');
            dbf = file.substr(0, dot_index);
        }
        dbf += ".dbf";
        dbfHandle = DBFCreate(dbf.c_str());
    }

    //create nodenumber field
    int node_field_id = -1;
    if (!isTxt) node_field_id = DBFAddField(dbfHandle, "nodenumber", FTInteger, 12, 0);
    vector<int> intFldIDs; //create integer fields
    vector<int> dbFldIDs;  //create double fields
    vector<int> strFldIDs; //create string fields
    vector<int> logFldIDs; //log fields
    if (!isTxt)
    {
        addDBHandleFields(dbfHandle, intFldNum, intFldNames, doubleFldNum, dbFldNames, intFldIDs,
                          dbFldIDs);
        addDBHandleStrFields(dbfHandle, strFldNum, strFldNames, strFldIDs);
        addDBHandleStrFields(dbfHandle, logFldNum, logFldNames, logFldIDs);
    }
    else
    {
        ofile << "nodenumber" << _IFO_DEL;
        for (int i = 0; i < intFldNum; ++i)
        {
            ofile << intFldNames[i] << _IFO_DEL;
        }
        for (int i = 0; i < doubleFldNum; ++i)
        {
            ofile << dbFldNames[i] << _IFO_DEL;
        }
        for (int i = 0; i < strFldNum; i++)
        {
            ofile << strFldNames[i].c_str() << _IFO_DEL;
        }
        for (int i = 0; i < logFldNum; i++)
        {
            ofile << logFldNames[i].c_str() << _IFO_DEL;
        }
        ofile << endl;
    }

    vector<vector<int>>::iterator iit = intfldVals.begin();
    vector<vector<double>>::iterator dit = dbfldVals.begin();
    vector<vector<string>>::iterator sit = strfldVals.begin();
    vector<vector<string>>::iterator lit = logfldVals.begin();
    vector<double> xs, ys; //vertex coordinates of the shape being written
    double *x = nullptr, *y = nullptr;
    //write points
    if (wtype == WRITE_POINT)
    {
        xs.assign(1, 0.0);
        ys.assign(1, 0.0);
        x = xs.data();
        y = ys.data();
        vector<PointCellIntersectionInfo> & pntInters =
            grid_intersection->getPointIntersectionInfo();
        vector<PointCellIntersectionInfo>::iterator pit = pntInters.begin();
        for (; pit != pntInters.end(); ++pit, ++iit, ++dit, ++sit, ++lit)
        {
            vector<int> & intVals = *iit;
            vector<double> & dbVals = *dit;
            vector<string> & strVals = *sit;
            vector<string> & logVals = *lit;
            x[0] = pit->point[0];
            y[0] = pit->point[1];

            SHPObject * shpObj = nullptr;
            int entity = -1;
            if (!isTxt)
            {
                shpObj = SHPCreateSimpleObject(SHPT_POINT, 1, x, y, nullptr);
                entity = SHPWriteObject(handler, -1, shpObj);
            }

            //write nodenumber
            if (!isTxt)
                DBFWriteIntegerAttribute(dbfHandle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1) << _IFO_DEL;
            writeIntDBVal(dbfHandle, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs,
                          strVals, logFldIDs, logVals, ofile, isTxt);
            if (!isTxt) SHPDestroyObject(shpObj);
        }
    }
    else if (wtype == WRITE_POLYLINE)
    {
        vector<PolylineCellIntersectionInfo> & plineInters =
            grid_intersection->getPolylineIntersectionInfo();
        vector<PolylineCellIntersectionInfo>::iterator pit = plineInters.begin();
        for (; pit != plineInters.end(); ++pit, ++iit, ++dit, ++sit, ++lit)
        {
            PolylineCellIntersectionInfo & lineInfo = *pit;
            vector<int> & intVals = *iit;
            vector<double> & dbVals = *dit;
            vector<string> & strVals = *sit;
            vector<string> & logVals = *lit;

            int vsize = lineInfo.m_plyline->getSize();
            if (vsize == 1)
            {
                cerr << "! Warning: bad polyline, ignored!" << endl;
                continue;
            }
            xs.assign(vsize, 0.0);
            ys.assign(vsize, 0.0);
            x = xs.data();
            y = ys.data();
            ply_vertex * pv = lineInfo.m_plyline->getHead();
            x[0] = pv->getPos()[0];
            y[0] = pv->getPos()[1];
            int vertex_index = 0;
            do
            {
                pv = pv->getNext();
                vertex_index++;
                x[vertex_index] = pv->getPos()[0];
                y[vertex_index] = pv->getPos()[1];
            } while (pv != lineInfo.m_plyline->getTail());

            SHPObject * shpObj = nullptr;
            int entity = -1;
            if (!isTxt)
            {
                shpObj = SHPCreateSimpleObject(SHPT_ARC, vsize, x, y, nullptr);
                entity = SHPWriteObject(handler, -1, shpObj);
            }

            //write nodenumber
            if (!isTxt)
                DBFWriteIntegerAttribute(dbfHandle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1) << _IFO_DEL;
            writeIntDBVal(dbfHandle, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs,
                          strVals, logFldIDs, logVals, ofile, isTxt);

            if (!isTxt) SHPDestroyObject(shpObj);
        }
    }
    else if (wtype == WRITE_POLYGON)
    {
        vector<PolygonCellIntersectionInfo> & pgonInters =
            grid_intersection->getPolygonIntersectionInfo();
        vector<PolygonCellIntersectionInfo>::iterator pit = pgonInters.begin();
        for (int iShape = 0; pit != pgonInters.end(); ++pit, ++iit, ++dit, ++iShape, ++sit, ++lit)
        {
            PolygonCellIntersectionInfo & pgonInfo = *pit;
            vector<int> & intVals = *iit;
            vector<double> & dbVals = *dit;
            vector<string> & strVals = *sit;
            vector<string> & logVals = *lit;

            int partSize = pgonInfo.m_interPolygon->size();
            int vsize = pgonInfo.m_interPolygon->getSize();
            if (vsize < 3)
            {
                cerr << "! Warning: bad polygon, ignored!" << endl;
                continue;
            }
            // each ring repeats its first vertex at the end
            vsize += partSize;
            xs.assign(vsize, 0.0);
            ys.assign(vsize, 0.0);
            x = xs.data();
            y = ys.data();

            //make sure the out boundary
            list<c_polygon::iterator> orderedPItrs;
            for (c_polygon::iterator plyit = pgonInfo.m_interPolygon->begin();
                 plyit != pgonInfo.m_interPolygon->end(); ++plyit)
            {
                if (plyit->getType() == c_ply::POUT)
                    orderedPItrs.push_front(plyit);
                else
                    orderedPItrs.push_back(plyit);
            }

            //part number
            int part_start = 0;
            vector<int> parts(partSize, 0);
            int * nParts = parts.data();
            int k = 0, vid = 0;
            for (list<c_polygon::iterator>::iterator lit = orderedPItrs.begin();
                 lit != orderedPItrs.end(); ++lit, k++)
            {
                c_polygon::iterator & plyit = *lit;
                nParts[k] = part_start;
                part_start += plyit->getSize() + 1;

                // the shapefile specification orders the outer ring clockwise
                // and holes counterclockwise
                ply_vertex * head = plyit->getHead();
                bool forward = (ringArea(head) < 0.0) == (plyit->getType() == c_ply::POUT);
                ply_vertex * pv = head;
                do
                {
                    x[vid] = pv->getPos()[0];
                    y[vid] = pv->getPos()[1];
                    pv = forward ? pv->getNext() : pv->getPre();
                    vid++;
                } while (pv != head);
                x[vid] = head->getPos()[0];
                y[vid] = head->getPos()[1];
                vid++;
            }

            //create a polygon shape here
            //SHPObject * shpObj = SHPCreateObject(SHPT_POLYGON,vsize,x, y,NULL);
            SHPObject * shpObj = nullptr;
            int entity = 0;
            if (!isTxt)
            {
                shpObj = SHPCreateObject(SHPT_POLYGON, iShape, partSize, nParts, nullptr, vsize, x,
                                         y, nullptr, nullptr);
                entity = SHPWriteObject(handler, -1, shpObj);
            }
            //write nodenumber
            if (!isTxt)
                DBFWriteIntegerAttribute(dbfHandle, entity, node_field_id,
                                         pit->box->number == -1 ? -1 : pit->box->number + 1);
            else
                ofile << (pit->box->number == -1 ? -1 : pit->box->number + 1) << _IFO_DEL;
            writeIntDBVal(dbfHandle, entity, intFldIDs, intVals, dbFldIDs, dbVals, strFldIDs,
                          strVals, logFldIDs, logVals, ofile, isTxt);

            if (!isTxt) SHPDestroyObject(shpObj);
        }
    }

    if (!isTxt)
    {
        DBFClose(dbfHandle);
        SHPClose(handler);
    }
    else
        ofile.close();
}

//write the grid and point intersection
void writePntIntersect(GridIntersection * grid_intersection, string file,
                       vector<string> & attributes, FILETYPE filetype)
{
    vector<PointCellIntersectionInfo> & pntInters = grid_intersection->getPointIntersectionInfo();

    vector<string> intFldName, dbFldName, strFldName, logFldName;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

    //parse the attributes
    bool write_point_id =
        true; // checkAttributesContain("pointid", attributes);//write shapefile_pointid
    if (write_point_id) intFldName.push_back(string("pointid"));
    //add the attributes in shape file
    checkShpAttrContain(grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld, intFldName, dbFldName,
                        strFldName, logFldName, attributes);

    for (auto & info : pntInters)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;
        if (write_point_id) int_values.push_back(info.m_pointID);

        addShpAttrVal(info.m_pointID, grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld,
                      int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool isTxt = (filetype == WRITE_TXT);
        writeIntersection(file, intFldName.size(), intflds, intFldName, dbFldName.size(), dbflds,
                          dbFldName, strFldName.size(), strflds, strFldName, logFldName.size(),
                          logflds, logFldName, grid_intersection, WRITE_POINT, isTxt);
    }
    else
    {
    }
}

//write the line and grid intersection
void writeLineIntersect(GridIntersection * grid_intersection, string file,
                        vector<string> & attributes, FILETYPE filetype)
{
    vector<PolylineCellIntersectionInfo> & polyline_infos =
        grid_intersection->getPolylineIntersectionInfo();

    vector<string> intFldName, dbFldName, strFldName, logFldName;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

    bool write_arc_id = true; // checkAttributesContain("arcid", attributes);//write shapefile_arcid
    bool write_length = true; //checkAttributesContain("length", attributes);//write length
    bool write_start_distance =
        true; //checkAttributesContain("starting_distance", attributes);//write starting_distance
    bool write_end_distance =
        true; //checkAttributesContain("ending_distance", attributes);//write ending_distance

    if (write_arc_id) intFldName.push_back(string("arcid"));
    if (write_length) dbFldName.push_back(string("length"));
    if (write_start_distance) dbFldName.push_back(string("starting_distance"));
    if (write_end_distance) dbFldName.push_back(string("ending_distance"));
    //add the shape file attributes
    checkShpAttrContain(grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld, intFldName, dbFldName,
                        strFldName, logFldName, attributes);

    for (auto & info : polyline_infos)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;

        if (write_arc_id) int_values.push_back(info.m_curveID);
        if (write_length) dbs.push_back(info.m_len);
        if (write_start_distance) dbs.push_back(info.m_startGeoDist);
        if (write_end_distance) dbs.push_back(info.m_endGeoDist);

        addShpAttrVal(info.m_curveID, grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld,
                      int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool isTxt = (filetype == WRITE_TXT);
        writeIntersection(file, intFldName.size(), intflds, intFldName, dbFldName.size(), dbflds,
                          dbFldName, strFldName.size(), strflds, strFldName, logFldName.size(),
                          logflds, logFldName, grid_intersection, WRITE_POLYLINE, isTxt);
    }
    else
    {
    }
}

//writhe the polygon and grid intersection
void writePolyIntersect(GridIntersection * grid_intersection, string file,
                        vector<string> & attributes, FILETYPE filetype)
{
    vector<PolygonCellIntersectionInfo> & sci = grid_intersection->getPolygonIntersectionInfo();

    vector<string> intFldName, dbFldName, strFldName, logFldName;
    vector<vector<int>> intflds;
    vector<vector<double>> dbflds;
    vector<vector<string>> strflds;
    vector<vector<string>> logflds;
    //to record whether each attribute needs to be written
    vector<bool> bIntFld, bDBFld, bStrFld, bLogFld;

    bool write_polygon_id =
        true;               // checkAttributesContain("polyid", attributes);//write shapefile_polyid
    bool write_area = true; //checkAttributesContain("totalarea", attributes);//write totalarea
    if (write_polygon_id) intFldName.push_back(string("polyid"));
    if (write_area) dbFldName.push_back(string("totalarea"));

    checkShpAttrContain(grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld, intFldName, dbFldName,
                        strFldName, logFldName, attributes);

    for (auto & info : sci)
    {
        vector<int> int_values;
        vector<double> dbs;
        vector<string> strs;
        vector<string> logs;

        if (write_polygon_id) int_values.push_back(info.m_polyID);
        if (write_area) dbs.push_back(info.m_area);

        addShpAttrVal(info.m_polyID, grid_intersection, bIntFld, bDBFld, bStrFld, bLogFld,
                      int_values, dbs, strs, logs);

        intflds.push_back(int_values);
        dbflds.push_back(dbs);
        strflds.push_back(strs);
        logflds.push_back(logs);
    }

    if (filetype == WRITE_SHP || filetype == WRITE_TXT)
    {
        bool isTxt = (filetype == WRITE_TXT);
        writeIntersection(file, intFldName.size(), intflds, intFldName, dbFldName.size(), dbflds,
                          dbFldName, strFldName.size(), strflds, strFldName, logFldName.size(),
                          logflds, logFldName, grid_intersection, WRITE_POLYGON, isTxt);
    }
    else
    {
    }
}
